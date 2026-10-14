#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Required package: data.table")
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Required package: fgsea")
})

library(data.table)
library(fgsea)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

DUR <- file.path(
  ROOT, "05_meta/moderators_v1_2",
  "HFD_duration_meta_regression_v1_2.tsv.gz"
)

MULTI <- file.path(
  ROOT, "05_meta/moderator_followup_v1",
  "HFD_duration_plus_fat_multivariable_meta_regression_v1.tsv.gz"
)

ANNOT <- file.path(
  ROOT, "05_meta/gene_level_v1_2",
  "HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz"
)

HALLMARK <- file.path(
  ROOT, "00_metadata/analysis_inputs/pathway_analysis_v1_2",
  "MSigDB_mouse_Hallmark_frozen_v1_2.tsv.gz"
)

OUTDIR <- file.path(ROOT, "06_pathway", "duration_modifier_v1")
METADIR <- file.path(
  ROOT, "00_metadata/analysis_outputs/duration_pathway_v1"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT_RANK_PRIMARY <- file.path(
  OUTDIR, "HFD_duration_primary_ranked_genes_z_v1.tsv.gz"
)
OUT_RANK_MULTI <- file.path(
  OUTDIR, "HFD_duration_adjusted_for_fat_ranked_genes_z_v1.tsv.gz"
)
OUT_GSEA_PRIMARY <- file.path(
  OUTDIR, "HFD_duration_Hallmark_GSEA_primary_v1.tsv"
)
OUT_GSEA_MULTI <- file.path(
  OUTDIR, "HFD_duration_Hallmark_GSEA_adjusted_for_fat_v1.tsv"
)
OUT_COMPARE <- file.path(
  OUTDIR, "HFD_duration_Hallmark_primary_vs_adjusted_comparison_v1.tsv"
)
OUT_LE <- file.path(
  OUTDIR, "HFD_duration_Hallmark_primary_leading_edge_v1.tsv.gz"
)
OUT_SUMMARY <- file.path(
  METADIR, "HFD_duration_Hallmark_GSEA_summary_v1.txt"
)

for (q in c(DUR, MULTI, ANNOT, HALLMARK)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

dur <- fread(DUR)
multi <- fread(MULTI)
ann <- fread(ANNOT)
hall <- fread(HALLMARK)

need_dur <- c("Geneid", "beta", "SE", "z_value", "p_value", "FDR_BH",
              "N_independent_GSE")
need_multi <- c("Geneid", "beta_duration_per4", "SE_duration_per4",
                "z_duration_per4", "p_duration_per4", "duration_FDR_BH",
                "N_independent_GSE")
missing <- setdiff(need_dur, names(dur))
if (length(missing)) stop("Duration table missing: ", paste(missing, collapse=", "))
missing <- setdiff(need_multi, names(multi))
if (length(missing)) stop("Multivariable table missing: ", paste(missing, collapse=", "))

ann2 <- unique(ann[, .(Geneid, gene_symbol, gene_type)])
if (anyDuplicated(ann2$Geneid)) stop("Duplicate Geneid in annotation table")

dur <- merge(dur, ann2, by="Geneid", all.x=TRUE, sort=FALSE)
multi <- merge(multi, ann2, by="Geneid", all.x=TRUE, sort=FALSE)

dur[, ensembl_gene := sub("\\.[0-9]+$", "", Geneid)]
multi[, ensembl_gene := sub("\\.[0-9]+$", "", Geneid)]

if (anyDuplicated(dur$ensembl_gene)) stop("Duplicate stripped Ensembl ID in duration table")
if (anyDuplicated(multi$ensembl_gene)) stop("Duplicate stripped Ensembl ID in multivariable table")

dur[, duration_rank_z := as.numeric(z_value)]
multi[, duration_rank_z_adjusted := as.numeric(z_duration_per4)]

if (any(!is.finite(dur$duration_rank_z))) stop("Non-finite duration z statistic")
if (any(!is.finite(multi$duration_rank_z_adjusted))) stop("Non-finite adjusted duration z statistic")

dur_rank <- dur[, .(
  Geneid, ensembl_gene, gene_symbol, gene_type,
  duration_rank_z,
  beta_duration_per4weeks = beta,
  SE_duration_per4weeks = SE,
  p_value,
  FDR_BH,
  N_independent_GSE
)]
setorder(dur_rank, -duration_rank_z, ensembl_gene)
fwrite(dur_rank, OUT_RANK_PRIMARY, sep="\t", quote=FALSE, na="NA", compress="gzip")

multi_rank <- multi[, .(
  Geneid, ensembl_gene, gene_symbol, gene_type,
  duration_rank_z_adjusted,
  beta_duration_per4weeks_adjusted_for_fat = beta_duration_per4,
  SE_duration_per4weeks_adjusted_for_fat = SE_duration_per4,
  p_duration_per4,
  duration_FDR_BH,
  N_independent_GSE
)]
setorder(multi_rank, -duration_rank_z_adjusted, ensembl_gene)
fwrite(multi_rank, OUT_RANK_MULTI, sep="\t", quote=FALSE, na="NA", compress="gzip")

pathways_all <- split(hall$ensembl_gene, hall$gs_name)

run_fgsea <- function(rank_dt, stat_col, label) {
  stats <- rank_dt[[stat_col]]
  names(stats) <- rank_dt$ensembl_gene
  stats <- sort(stats, decreasing=TRUE)

  overlaps <- sapply(pathways_all, function(g) sum(unique(g) %in% names(stats)))
  keep <- names(overlaps)[overlaps >= 15 & overlaps <= 500]
  pathways <- pathways_all[keep]

  if (length(pathways) == 0L) stop(label, ": no eligible Hallmark pathways")

  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=stats,
    minSize=15,
    maxSize=500,
    eps=0,
    scoreType="std",
    nproc=1,
    gseaParam=1
  )

  r <- as.data.table(fg)
  r[, padj_BH := p.adjust(pval, method="BH")]
  r[, leadingEdge_ensembl := vapply(
    leadingEdge, function(z) paste(z, collapse=";"), character(1)
  )]
  r[, leadingEdge_size := lengths(leadingEdge)]

  r[, abs_NES_sort := abs(NES)]
  setorder(r, padj_BH, -abs_NES_sort, pathway)
  r[, abs_NES_sort := NULL]

  r[, Analysis := label]
  r[]
}

set.seed(20260904)
g_primary <- run_fgsea(dur_rank, "duration_rank_z", "PRIMARY_DURATION")

set.seed(20260904)
g_multi <- run_fgsea(
  multi_rank,
  "duration_rank_z_adjusted",
  "DURATION_ADJUSTED_FOR_FAT"
)

fwrite(
  g_primary[, !"leadingEdge", with=FALSE],
  OUT_GSEA_PRIMARY,
  sep="\t", quote=FALSE, na="NA"
)
fwrite(
  g_multi[, !"leadingEdge", with=FALSE],
  OUT_GSEA_MULTI,
  sep="\t", quote=FALSE, na="NA"
)

cmp <- merge(
  g_primary[, .(
    pathway,
    NES_primary=NES,
    p_primary=pval,
    FDR_primary=padj_BH
  )],
  g_multi[, .(
    pathway,
    NES_adjusted_for_fat=NES,
    p_adjusted_for_fat=pval,
    FDR_adjusted_for_fat=padj_BH
  )],
  by="pathway",
  all=TRUE
)
cmp[, Same_NES_direction := ifelse(
  sign(NES_primary) == sign(NES_adjusted_for_fat), "YES", "NO"
)]
cmp[, abs_NES_primary_sort := abs(NES_primary)]
setorder(cmp, FDR_primary, -abs_NES_primary_sort, pathway)
cmp[, abs_NES_primary_sort := NULL]
fwrite(cmp, OUT_COMPARE, sep="\t", quote=FALSE, na="NA")

sig_paths <- g_primary[padj_BH < 0.05, pathway]
if (length(sig_paths)) {
  le <- g_primary[pathway %in% sig_paths, .(
    ensembl_gene=unlist(leadingEdge)
  ), by=.(pathway, NES, pval, padj=padj_BH)]

  le <- merge(
    le,
    dur_rank[, .(
      ensembl_gene, Geneid, gene_symbol, gene_type,
      duration_rank_z,
      beta_duration_per4weeks,
      FDR_BH,
      N_independent_GSE
    )],
    by="ensembl_gene",
    all.x=TRUE,
    sort=FALSE
  )

  le[, abs_rank_sort := abs(duration_rank_z)]
  setorder(le, pathway, -abs_rank_sort, ensembl_gene)
  le[, abs_rank_sort := NULL]
} else {
  le <- data.table(
    pathway=character(), NES=numeric(), pval=numeric(), padj=numeric(),
    ensembl_gene=character(), Geneid=character(), gene_symbol=character(),
    gene_type=character(), duration_rank_z=numeric(),
    beta_duration_per4weeks=numeric(), FDR_BH=numeric(),
    N_independent_GSE=integer()
  )
}
fwrite(le, OUT_LE, sep="\t", quote=FALSE, na="NA", compress="gzip")

n_sig_primary <- sum(g_primary$padj_BH < .05, na.rm=TRUE)
n_pos_primary <- sum(g_primary$padj_BH < .05 & g_primary$NES > 0, na.rm=TRUE)
n_neg_primary <- sum(g_primary$padj_BH < .05 & g_primary$NES < 0, na.rm=TRUE)

n_sig_multi <- sum(g_multi$padj_BH < .05, na.rm=TRUE)
n_pos_multi <- sum(g_multi$padj_BH < .05 & g_multi$NES > 0, na.rm=TRUE)
n_neg_multi <- sum(g_multi$padj_BH < .05 & g_multi$NES < 0, na.rm=TRUE)

common <- cmp[is.finite(NES_primary) & is.finite(NES_adjusted_for_fat)]
pearson_nes <- cor(common$NES_primary, common$NES_adjusted_for_fat, method="pearson")
spearman_nes <- cor(common$NES_primary, common$NES_adjusted_for_fat, method="spearman")
same_dir <- mean(common$Same_NES_direction == "YES")

top_pos <- g_primary[NES > 0][order(padj_BH, -NES)]
top_neg <- g_primary[NES < 0][order(padj_BH, NES)]

fmt_top <- function(x, n=10L) {
  if (!nrow(x)) return("NONE")
  x <- head(x, n)
  paste(
    paste0(
      x$pathway,
      "|NES=", format(round(x$NES, 3), nsmall=3),
      "|FDR=", format(x$padj_BH, digits=3, scientific=TRUE)
    ),
    collapse="\n"
  )
}

summary <- c(
  "===== HFD DURATION HALLMARK GSEA SUMMARY v1 =====",
  "",
  "STATUS=PASS",
  paste0("R_version=", R.version.string),
  paste0("fgsea_version=", as.character(packageVersion("fgsea"))),
  "MSigDB_release=2025.1.Mm",
  "MSigDB_collection=MH_Hallmark",
  "",
  "PRIMARY_RANK=duration_meta_regression_z_statistic",
  "PRIMARY_EFFECT_UNIT=change_in_HFD_log2FC_per_4_weeks",
  paste0("Primary_ranked_genes=", nrow(dur_rank)),
  paste0("Primary_Hallmark_pathways_tested=", nrow(g_primary)),
  paste0("Primary_Hallmark_FDR_lt_0.05=", n_sig_primary),
  paste0("Primary_positive_NES_FDR05=", n_pos_primary),
  paste0("Primary_negative_NES_FDR05=", n_neg_primary),
  "",
  "SENSITIVITY_RANK=duration_z_adjusted_for_fat_percent",
  paste0("Sensitivity_ranked_genes=", nrow(multi_rank)),
  paste0("Sensitivity_Hallmark_pathways_tested=", nrow(g_multi)),
  paste0("Sensitivity_Hallmark_FDR_lt_0.05=", n_sig_multi),
  paste0("Sensitivity_positive_NES_FDR05=", n_pos_multi),
  paste0("Sensitivity_negative_NES_FDR05=", n_neg_multi),
  "",
  paste0("Primary_vs_adjusted_NES_Pearson_r=", format(pearson_nes, digits=6)),
  paste0("Primary_vs_adjusted_NES_Spearman_rho=", format(spearman_nes, digits=6)),
  paste0("Primary_vs_adjusted_same_direction_fraction=", format(same_dir, digits=6)),
  "",
  "TOP_PRIMARY_POSITIVE_DURATION_PATHWAYS",
  fmt_top(top_pos, 10L),
  "",
  "TOP_PRIMARY_NEGATIVE_DURATION_PATHWAYS",
  fmt_top(top_neg, 10L),
  "",
  paste0("Primary_rank=", sub(paste0("^",ROOT,"/"),"",OUT_RANK_PRIMARY)),
  paste0("Adjusted_rank=", sub(paste0("^",ROOT,"/"),"",OUT_RANK_MULTI)),
  paste0("Primary_GSEA=", sub(paste0("^",ROOT,"/"),"",OUT_GSEA_PRIMARY)),
  paste0("Adjusted_GSEA=", sub(paste0("^",ROOT,"/"),"",OUT_GSEA_MULTI)),
  paste0("Comparison=", sub(paste0("^",ROOT,"/"),"",OUT_COMPARE)),
  paste0("Primary_leading_edge=", sub(paste0("^",ROOT,"/"),"",OUT_LE)),
  "",
  "CAUSAL_INTERPRETATION=NO",
  "CONTROL_DIET_PATHWAY_ANALYSIS_PERFORMED=NO",
  "HETEROGENEITY_RANKED_PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, OUT_SUMMARY)
cat("\n", paste(summary, collapse="\n"), "\n", sep="")

