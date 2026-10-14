#!/usr/bin/env Rscript
# IJMS Reviewer 1, Major Comment 5: GSE-level aggregation sensitivity.
# Alternative analysis: one inverse-variance-weighted effect per contributing GSE,
# followed by a single-level REML random-effects meta-analysis per gene.
# IMPORTANT: its PI describes a new GSE-level aggregate, NOT a new individual
# contrast; PI-based gene membership is therefore not strictly equivalent
# to the primary two-level model's new-contrast PI.
suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
})
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[1], mustWork=TRUE) else getwd()
long <- file.path(root, '00_metadata/analysis_inputs/meta_analysis_v1',
                  'HFD_primary105_meta_effects_primaryQC_v1.tsv.gz')
primary <- file.path(root, '05_meta/gene_level_v1_2',
                     'HFD_gene_multilevel_meta_primary_v1_2.tsv.gz')
outdir <- file.path(root, '00_metadata/analysis_outputs/revision_R1C5_GSE_aggregate')
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
if (!file.exists(long) || !file.exists(primary)) {
  stop('Required primary input missing: ', long, '\n or ', primary)
}
cat('===== IJMS R1 Comment 5: GSE-aggregated alternative model =====\n')
cat('Project:', root, '\n')
cat('metafor:', as.character(packageVersion('metafor')), '\n')
# Full primary gene universe is retained for genome-wide BH adjustment.
p <- fread(primary, select=c('Geneid','pooled_log2FC','pooled_FDR_BH',
                              'stable_core_preLOO'))
if (anyDuplicated(p$Geneid)) stop('Primary results contain duplicate Geneid')
if (nrow(p)!=18234L) stop('Unexpected primary gene universe size: ', nrow(p))
p[, primary_core := toupper(as.character(stable_core_preLOO))=='YES']
if (sum(p$primary_core)!=59L) stop('Expected exactly 59 primary core genes')
cat('Reading contrast-level meta input...\n')
d <- fread(long, select=c('Geneid','GSE','Contrast_id','log2FoldChange',
                           'sampling_variance'))
d <- d[Geneid %chin% p$Geneid]
if (!nrow(d)) stop('No gene-level contrasts matched primary output')
if (any(!is.finite(d$log2FoldChange)) ||
    any(!is.finite(d$sampling_variance) | d$sampling_variance<=0)) {
  stop('Invalid contrast-level effect or sampling variance')
}
if (anyDuplicated(d[, .(Geneid, Contrast_id)])) {
  stop('Duplicate gene-contrast entries in primary input')
}
cat('Aggregating contrasts within each GSE...\n')
g <- d[, .(
  gse_yi = sum(log2FoldChange/sampling_variance)/sum(1/sampling_variance),
  gse_vi = 1/sum(1/sampling_variance),
  contrasts_in_GSE = .N
), by=.(Geneid,GSE)]
if (any(!is.finite(g$gse_yi) | !is.finite(g$gse_vi) | g$gse_vi<=0)) {
  stop('Invalid aggregated GSE-level effects')
}
setkey(g,Geneid)
groups <- split(g, by='Geneid',keep.by=FALSE)
cat('Fitting REML random-effects meta-analysis for ',length(groups),' genes...\n',sep='')
fit_one <- function(x, gene) {
  n <- nrow(x)
  if (n<10L) return(data.table(Geneid=gene,N_GSE=n,status='COVERAGE_FAIL'))
  fit <- tryCatch(
    suppressWarnings(rma.uni(yi=x$gse_yi,vi=x$gse_vi,
                             method='REML',test='z')),
    error=function(e)e)
  if (inherits(fit,'error'))
    return(data.table(Geneid=gene,N_GSE=n,status='FIT_FAIL',error=conditionMessage(fit)))
  mu <- as.numeric(coef(fit)[1]); se <- as.numeric(fit$se[1]);
  tau2 <- as.numeric(fit$tau2); pv <- as.numeric(fit$pval[1])
  l <- mu-qnorm(.975)*sqrt(se^2+tau2)
  u <- mu+qnorm(.975)*sqrt(se^2+tau2)
  con <- mean(sign(x$gse_yi)==sign(mu))
  data.table(Geneid=gene,N_GSE=n,status='PASS',
             alt_pooled_log2FC=mu,alt_SE=se,alt_p_value=pv,
             alt_tau2_between_GSE=tau2,alt_PI95_GSE_lower=l,
             alt_PI95_GSE_upper=u,alt_PI95_GSE_excludes_zero=(l>0 | u<0),
             alt_GSE_direction_consistency=con)
}
res <- rbindlist(Map(fit_one,groups,names(groups)), fill=TRUE)
if (anyDuplicated(res$Geneid)) stop('Duplicate genes in alternative results')
res[, alt_FDR_BH:=NA_real_]
res[status=='PASS',alt_FDR_BH:=p.adjust(alt_p_value,method='BH')]
# Missing model fits are excluded from BH denominator, so report any failures.
if (any(res$status!='PASS'))
  warning('Alternative fits failed; examine failures and adjust FDR universe interpretation')
setkey(res,Geneid);setkey(p,Geneid)
joined <- merge(p,res,by='Geneid',all.x=TRUE,sort=FALSE)
joined[, alt_same_pooled_sign := status=='PASS' &
       sign(pooled_log2FC)==sign(alt_pooled_log2FC)]
joined[, alt_FDR_sig := status=='PASS' & alt_FDR_BH<.05]
joined[, alt_exploratory_core := alt_FDR_sig &
       alt_PI95_GSE_excludes_zero & alt_GSE_direction_consistency>=.80]
# Explicitly label exploratory alternative core: prediction targets differ.
file <- file.path(outdir,'GSE_aggregated_alternative_model_gene_results_v1.tsv.gz')
fwrite(joined,file,sep='\t')
num <- function(mask) sum(mask,na.rm=TRUE)
cs <- joined[primary_core==TRUE]
ss <- joined[status=='PASS' & is.finite(pooled_log2FC) &
               is.finite(alt_pooled_log2FC)]
summary <- c(
  '===== GSE-aggregated alternative REML sensitivity (R1C5) =====',
  paste('Original genes:',nrow(joined)),
  paste('Successful alternative fits:',num(joined$status=='PASS')),
  paste('Failed alternative fits:',num(joined$status!='PASS')),
  paste('Primary FDR<0.05:',num(joined$pooled_FDR_BH<.05)),
  paste('Alternative FDR<0.05:',num(joined$alt_FDR_sig)),
  paste('Original core genes:',nrow(cs)),
  paste('Original core with same sign under alternative:',num(cs$alt_same_pooled_sign)),
  paste('Original core FDR<0.05 under alternative:',num(cs$alt_FDR_sig)),
  paste('Original core with alternative GSE-level 95% PI excluding zero:',
        num(cs$alt_PI95_GSE_excludes_zero)),
  paste('Original core satisfying all three alternative criteria (exploratory):',
        num(cs$alt_exploratory_core)),
  paste('Entire exploratory alternative core:',num(joined$alt_exploratory_core)),
  paste('Primary versus alternative pooled LFC Pearson correlation:',
        round(cor(ss$pooled_log2FC,ss$alt_pooled_log2FC,method='pearson'),5)),
  paste('Primary versus alternative pooled LFC Spearman correlation:',
        round(cor(ss$pooled_log2FC,ss$alt_pooled_log2FC,method='spearman'),5)),
  'Note: The alternative 95% PI targets a new GSE-level *aggregated* effect.',
  'It is NOT directly equivalent to the primary model PI for a new contrast in a new GSE.',
  'Original core membership is never overwritten by this sensitivity analysis.'
)
writeLines(summary,file.path(outdir,'GSE_aggregated_sensitivity_summary_v1.txt'))
cat(paste(summary,collapse='\n'),'\n')
cat('Gene-level output: ',file,'\n',sep='')
