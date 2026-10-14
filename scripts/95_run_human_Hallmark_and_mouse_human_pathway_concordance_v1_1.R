#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table","fgsea","msigdbr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse=", "))
})
library(data.table)
library(fgsea)
library(msigdbr)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

HUMAN <- file.path(
  ROOT, "07_human_validation/gene_meta_v1",
  "human_MASLD_three_cohort_random_effects_meta_v1.tsv.gz"
)
HUMAN_STRICT <- file.path(
  ROOT, "07_human_validation/gene_meta_v1",
  "human_MASLD_three_cohort_strict4_sensitivity_meta_v1.tsv.gz"
)
MOUSE_PATH <- file.path(
  ROOT, "06_pathway/hallmark_mean_response_v1_2",
  "HFD_Hallmark_GSEA_pooled_meta_zstat_all_v1_2.tsv"
)

OUTDIR <- file.path(
  ROOT, "07_human_validation/pathway_concordance_v1_1"
)
METADIR <- file.path(
  ROOT,
  "00_metadata/human_validation/analysis_outputs/pathway_concordance_v1_1"
)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)

OUT_RE <- file.path(
  OUTDIR, "human_MASLD_Hallmark_GSEA_RE_zstat_v1_1.tsv"
)
OUT_FE <- file.path(
  OUTDIR, "human_MASLD_Hallmark_GSEA_FE_zstat_v1_1.tsv"
)
OUT_STRICT <- file.path(
  OUTDIR, "human_MASLD_Hallmark_GSEA_strict4_RE_zstat_v1_1.tsv"
)
OUT_JOIN <- file.path(
  OUTDIR, "mouse_HFD_human_MASLD_Hallmark_concordance_v1_1.tsv"
)
OUT_SUM <- file.path(
  METADIR, "mouse_human_Hallmark_concordance_summary_v1_1.txt"
)

for (q in c(HUMAN,HUMAN_STRICT,MOUSE_PATH)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

h <- fread(HUMAN)
hs <- fread(HUMAN_STRICT)
mp <- fread(MOUSE_PATH)

# ------------------------------------------------------------
# Human Hallmark collection, current local msigdbr package.
# Try current API first, then legacy API.
# ------------------------------------------------------------
hallmark <- tryCatch(
  msigdbr(species="Homo sapiens", collection="H"),
  error=function(e) NULL
)
if (is.null(hallmark)) {
  hallmark <- tryCatch(
    msigdbr(species="Homo sapiens", category="H"),
    error=function(e)e
  )
}
if (inherits(hallmark,"error"))
  stop("Could not load human Hallmark gene sets from msigdbr: ",
       conditionMessage(hallmark))

hallmark <- as.data.table(hallmark)

gs_col <- intersect(c("gs_name","gs_name"), names(hallmark))[1]
id_col <- intersect(
  c("ncbi_gene","entrez_gene","entrezgene_id","gene_id"),
  names(hallmark)
)[1]

if (is.na(gs_col) || is.na(id_col)) {
  stop(
    "Could not resolve Hallmark columns. Available: ",
    paste(names(hallmark), collapse=",")
  )
}

hallmark[, pathway := as.character(get(gs_col))]
hallmark[, ENTREZID := as.character(get(id_col))]
hallmark <- unique(hallmark[
  !is.na(pathway) & nzchar(pathway) &
  !is.na(ENTREZID) & nzchar(ENTREZID),
  .(pathway,ENTREZID)
])

pathways <- split(hallmark$ENTREZID, hallmark$pathway)

# ------------------------------------------------------------
# Human rankings.
# Primary = random-effects meta-analysis z statistic.
# Sensitivities = fixed-effect z, strict4 random-effects z.
# ------------------------------------------------------------
need_h <- c(
  "ENTREZID","human_RE_z","human_FE_z",
  "human_RE_log2FC","human_FE_log2FC"
)
miss <- setdiff(need_h,names(h))
if (length(miss))
  stop("Human primary meta missing columns: ",paste(miss,collapse=","))

need_hs <- c("ENTREZID","human_RE_z","human_RE_log2FC")
miss <- setdiff(need_hs,names(hs))
if (length(miss))
  stop("Human strict4 meta missing columns: ",paste(miss,collapse=","))

make_rank <- function(dt, idcol, statcol) {
  z <- dt[, .(
    ENTREZID=as.character(get(idcol)),
    stat=as.numeric(get(statcol))
  )]
  z <- z[is.finite(stat) & !is.na(ENTREZID) & nzchar(ENTREZID)]
  if (anyDuplicated(z$ENTREZID))
    stop("Duplicate ENTREZID in ranking: ", statcol)
  r <- z$stat
  names(r) <- z$ENTREZID
  sort(r,decreasing=TRUE)
}

rank_re <- make_rank(h,"ENTREZID","human_RE_z")
rank_fe <- make_rank(h,"ENTREZID","human_FE_z")
rank_strict <- make_rank(hs,"ENTREZID","human_RE_z")

run_fg <- function(rankvec, label) {
  set.seed(20260907)
  x <- fgseaMultilevel(
    pathways=pathways,
    stats=rankvec,
    minSize=15,
    maxSize=500,
    eps=0
  )
  x <- as.data.table(x)
  x[, Analysis := label]
  x[, pathway := as.character(pathway)]
  x[, leadingEdge := vapply(
    leadingEdge,
    function(v) paste(v,collapse=","),
    character(1)
  )]
  x[, abs_NES_sort := abs(NES)]
  setorder(x, padj, -abs_NES_sort, pathway)
  x[, abs_NES_sort := NULL]
  x
}

re <- run_fg(rank_re,"PRIMARY_human_RE_z")
fe <- run_fg(rank_fe,"SENSITIVITY_human_FE_z")
st <- run_fg(rank_strict,"SENSITIVITY_human_strict4_RE_z")

fwrite(re,OUT_RE,sep="\t",quote=FALSE,na="NA")
fwrite(fe,OUT_FE,sep="\t",quote=FALSE,na="NA")
fwrite(st,OUT_STRICT,sep="\t",quote=FALSE,na="NA")

# ------------------------------------------------------------
# Resolve mouse Hallmark columns robustly.
# ------------------------------------------------------------
resolve_one <- function(dt, exact, regex, label) {
  e <- exact[exact %in% names(dt)]
  if (length(e)) return(e[1])
  r <- grep(regex,names(dt),value=TRUE,ignore.case=TRUE,perl=TRUE)
  if (length(r)==1L) return(r)
  stop(
    "Could not uniquely resolve ",label,
    ". Candidates: ",paste(r,collapse=","),
    ". Available: ",paste(names(dt),collapse=",")
  )
}

m_path <- resolve_one(
  mp,
  c("pathway","Pathway","gs_name","term"),
  "(pathway|gs_name|term)",
  "mouse pathway"
)
m_nes <- resolve_one(
  mp,
  c("NES","nes"),
  "^NES$|^nes$",
  "mouse NES"
)
m_padj <- resolve_one(
  mp,
  c("padj","FDR","qval","qvalue"),
  "(padj|fdr|qval|qvalue)",
  "mouse adjusted P"
)

mouse <- mp[, .(
  pathway=as.character(get(m_path)),
  mouse_NES=as.numeric(get(m_nes)),
  mouse_FDR=as.numeric(get(m_padj))
)]

# Human compact tables.
hre <- re[, .(
  pathway,
  human_RE_NES=NES,
  human_RE_pval=pval,
  human_RE_FDR=padj,
  human_RE_size=size
)]
hfe <- fe[, .(
  pathway,
  human_FE_NES=NES,
  human_FE_FDR=padj
)]
hst <- st[, .(
  pathway,
  human_strict4_RE_NES=NES,
  human_strict4_RE_FDR=padj
)]

j <- Reduce(
  function(x,y) merge(x,y,by="pathway",all=FALSE),
  list(mouse,hre,hfe,hst)
)

if (nrow(j) < 45L)
  stop("Unexpectedly few Hallmark pathways in mouse-human join: ",nrow(j))

j[, same_direction_RE := sign(mouse_NES)==sign(human_RE_NES)]
j[, same_direction_FE := sign(mouse_NES)==sign(human_FE_NES)]
j[, same_direction_strict4 := sign(mouse_NES)==sign(human_strict4_RE_NES)]
j[, both_FDR05_RE := mouse_FDR < .05 & human_RE_FDR < .05]
j[, both_FDR05_same_direction_RE := both_FDR05_RE & same_direction_RE]

j[, abs_human_RE_sort := abs(human_RE_NES)]
setorder(j,human_RE_FDR,-abs_human_RE_sort)
j[, abs_human_RE_sort := NULL]

fwrite(j,OUT_JOIN,sep="\t",quote=FALSE,na="NA")

# ------------------------------------------------------------
# Concordance summaries.
# ------------------------------------------------------------
pear_re <- cor(j$mouse_NES,j$human_RE_NES,method="pearson")
spear_re <- cor(j$mouse_NES,j$human_RE_NES,method="spearman")
pear_fe <- cor(j$mouse_NES,j$human_FE_NES,method="pearson")
spear_fe <- cor(j$mouse_NES,j$human_FE_NES,method="spearman")
pear_st <- cor(j$mouse_NES,j$human_strict4_RE_NES,method="pearson")
spear_st <- cor(j$mouse_NES,j$human_strict4_RE_NES,method="spearman")

# Fisher z approximate CI for Pearson correlation across 50 pathways.
r_ci <- function(r,n) {
  z <- atanh(r)
  se <- 1/sqrt(n-3)
  c(low=tanh(z-1.96*se), high=tanh(z+1.96*se))
}
ci_re <- r_ci(pear_re,nrow(j))

summary <- c(
  "===== MOUSE HFD - HUMAN MASLD HALLMARK CONCORDANCE SUMMARY v1.1 =====",
  "",
  "STATUS=PASS",
  "SUPERSEDES=v1_data.table_setorder_abs_expression_bug",
  "ROLE=SECONDARY_CROSS_SPECIES_PATHWAY_VALIDATION",
  "PRIMARY_HUMAN_RANK=random_effects_meta_z_statistic",
  "SENSITIVITY_HUMAN_RANKS=fixed_effect_z_and_GSE130970_strict4_RE_z",
  "MOUSE_PATHWAY_INPUT=HFD_Hallmark_GSEA_pooled_meta_zstat_all_v1_2",
  "HALLMARK_SOURCE=msigdbr_Homo_sapiens_H",
  "GSEA_METHOD=fgseaMultilevel",
  "GSEA_MIN_SIZE=15",
  "GSEA_MAX_SIZE=500",
  "",
  paste0("Human_RE_ranked_genes=",length(rank_re)),
  paste0("Human_FE_ranked_genes=",length(rank_fe)),
  paste0("Human_strict4_ranked_genes=",length(rank_strict)),
  paste0("Hallmark_pathways_joined=",nrow(j)),
  "",
  paste0("Human_RE_Hallmark_FDR_lt_0.05=",sum(re$padj<.05,na.rm=TRUE)),
  paste0("Human_RE_Hallmark_FDR05_positive=",sum(re$padj<.05 & re$NES>0,na.rm=TRUE)),
  paste0("Human_RE_Hallmark_FDR05_negative=",sum(re$padj<.05 & re$NES<0,na.rm=TRUE)),
  "",
  paste0("Mouse_vs_human_RE_NES_Pearson_r=",format(pear_re,digits=6)),
  paste0(
    "Mouse_vs_human_RE_NES_Pearson_95CI=",
    format(ci_re["low"],digits=6),"_to_",format(ci_re["high"],digits=6)
  ),
  paste0("Mouse_vs_human_RE_NES_Spearman_rho=",format(spear_re,digits=6)),
  paste0("Mouse_vs_human_RE_same_direction_fraction=",format(mean(j$same_direction_RE),digits=6)),
  paste0("Mouse_and_human_RE_both_FDR05=",sum(j$both_FDR05_RE,na.rm=TRUE)),
  paste0(
    "Mouse_and_human_RE_both_FDR05_same_direction=",
    sum(j$both_FDR05_same_direction_RE,na.rm=TRUE)
  ),
  "",
  paste0("Mouse_vs_human_FE_NES_Pearson_r=",format(pear_fe,digits=6)),
  paste0("Mouse_vs_human_FE_NES_Spearman_rho=",format(spear_fe,digits=6)),
  paste0("Mouse_vs_human_FE_same_direction_fraction=",format(mean(j$same_direction_FE),digits=6)),
  "",
  paste0("Mouse_vs_human_strict4_RE_NES_Pearson_r=",format(pear_st,digits=6)),
  paste0("Mouse_vs_human_strict4_RE_NES_Spearman_rho=",format(spear_st,digits=6)),
  paste0(
    "Mouse_vs_human_strict4_RE_same_direction_fraction=",
    format(mean(j$same_direction_strict4),digits=6)
  ),
  "",
  paste0("Human_RE_GSEA=",sub(paste0("^",ROOT,"/"),"",OUT_RE)),
  paste0("Human_FE_GSEA=",sub(paste0("^",ROOT,"/"),"",OUT_FE)),
  paste0("Human_strict4_GSEA=",sub(paste0("^",ROOT,"/"),"",OUT_STRICT)),
  paste0("Mouse_human_join=",sub(paste0("^",ROOT,"/"),"",OUT_JOIN)),
  "",
  "CAUSAL_INTERPRETATION=NO",
  "GENE_LEVEL_PRIMARY_TRANSLATIONAL_INFERENCE_CHANGED=NO"
)

writeLines(summary,OUT_SUM)
cat("\n",paste(summary,collapse="\n"),"\n",sep="")

