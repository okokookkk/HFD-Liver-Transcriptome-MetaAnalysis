#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Required package: data.table")
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Required package: fgsea")
})
library(data.table)
library(fgsea)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

META <- file.path(
  ROOT, "05_meta/gene_level_v1_2",
  "HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz"
)
HALLMARK <- file.path(
  ROOT, "00_metadata/analysis_inputs/pathway_analysis_v1_2",
  "MSigDB_mouse_Hallmark_frozen_v1_2.tsv.gz"
)
MEAN_GSEA <- file.path(
  ROOT, "06_pathway/hallmark_mean_response_v1_2",
  "HFD_Hallmark_GSEA_pooled_meta_zstat_all_v1_2.tsv"
)
DURATION_GSEA <- file.path(
  ROOT, "06_pathway/duration_modifier_v1",
  "HFD_duration_Hallmark_GSEA_primary_v1.tsv"
)

OUTDIR <- file.path(ROOT, "06_pathway", "heterogeneity_v1")
METADIR <- file.path(
  ROOT, "00_metadata/analysis_outputs/heterogeneity_pathway_v1"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT_RANK_Q <- file.path(
  OUTDIR, "HFD_gene_heterogeneity_standardized_Q_rank_v1.tsv.gz"
)
OUT_Q <- file.path(
  OUTDIR, "HFD_Hallmark_heterogeneity_standardized_Q_GSEA_v1.tsv"
)
OUT_I2 <- file.path(
  OUTDIR, "HFD_Hallmark_heterogeneity_I2_GSEA_v1.tsv"
)
OUT_Q30 <- file.path(
  OUTDIR, "HFD_Hallmark_heterogeneity_standardized_Q_ge30GSE_GSEA_v1.tsv"
)
OUT_MAP <- file.path(
  OUTDIR, "HFD_Hallmark_mean_duration_heterogeneity_map_v1.tsv"
)
OUT_SUMMARY <- file.path(
  METADIR, "HFD_heterogeneity_Hallmark_GSEA_summary_v1.txt"
)

for (q in c(META, HALLMARK, MEAN_GSEA, DURATION_GSEA)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

meta <- fread(META)
hall <- fread(HALLMARK)
mean_gsea <- fread(MEAN_GSEA)
dur_gsea <- fread(DURATION_GSEA)

need_meta <- c(
  "Geneid", "gene_symbol", "gene_type",
  "N_independent_GSE", "QE", "QE_df", "QE_p_value",
  "tau2_total", "I2_total_percent"
)
missing <- setdiff(need_meta, names(meta))
if (length(missing))
  stop("Meta table missing columns: ", paste(missing, collapse=", "))

for (cc in c("N_independent_GSE","QE","QE_df","QE_p_value",
             "tau2_total","I2_total_percent")) {
  meta[, (cc) := as.numeric(get(cc))]
}

meta[, ensembl_gene := sub("\\.[0-9]+$", "", Geneid)]
if (anyDuplicated(meta$ensembl_gene))
  stop("Duplicate version-stripped Ensembl IDs")

if (any(!is.finite(meta$QE)) || any(!is.finite(meta$QE_df)) ||
    any(meta$QE_df <= 0))
  stop("Invalid QE/QE_df values")
if (any(!is.finite(meta$I2_total_percent)))
  stop("Non-finite I2_total_percent")

# Exploratory primary heterogeneity rank:
# standardized excess Cochran-type Q score.
# Positive = more residual heterogeneity than expected under homogeneity;
# negative = less. This standardizes for differing Q degrees of freedom.
meta[, heterogeneity_Q_z := (QE - QE_df) / sqrt(2 * QE_df)]

if (any(!is.finite(meta$heterogeneity_Q_z)))
  stop("Non-finite standardized heterogeneity Q score")

rank_q <- meta[, .(
  Geneid, ensembl_gene, gene_symbol, gene_type,
  N_independent_GSE,
  heterogeneity_Q_z,
  QE, QE_df, QE_p_value,
  tau2_total, I2_total_percent
)]
setorder(rank_q, -heterogeneity_Q_z, ensembl_gene)
fwrite(rank_q, OUT_RANK_Q, sep="\t", quote=FALSE, na="NA", compress="gzip")

pathways_all <- split(hall$ensembl_gene, hall$gs_name)

run_gsea <- function(ids, stats, score_type, label) {
  names(stats) <- ids
  ok <- is.finite(stats) & !is.na(names(stats)) & nzchar(names(stats))
  stats <- stats[ok]
  if (anyDuplicated(names(stats))) stop(label, ": duplicate named statistics")
  stats <- sort(stats, decreasing=TRUE)

  overlaps <- sapply(pathways_all, function(g) {
    sum(unique(g) %in% names(stats))
  })
  keep <- names(overlaps)[overlaps >= 15 & overlaps <= 500]
  pathways <- pathways_all[keep]

  if (!length(pathways)) stop(label, ": no eligible pathways")

  set.seed(20260904)
  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=stats,
    minSize=15,
    maxSize=500,
    eps=0,
    scoreType=score_type,
    nproc=1,
    gseaParam=1
  )

  r <- as.data.table(fg)
  r[, padj_BH := p.adjust(pval, method="BH")]
  r[, leadingEdge_ensembl := vapply(
    leadingEdge, function(z) paste(z, collapse=";"), character(1)
  )]
  r[, leadingEdge_size := lengths(leadingEdge)]
  r[, N_ranked_overlap := overlaps[pathway]]

  r[, abs_NES_sort := abs(NES)]
  setorder(r, padj_BH, -abs_NES_sort, pathway)
  r[, abs_NES_sort := NULL]
  r[, Analysis := label]
  r[]
}

# A) Primary exploratory: standardized excess-Q rank.
g_q <- run_gsea(
  meta$ensembl_gene,
  meta$heterogeneity_Q_z,
  "std",
  "STANDARDIZED_Q_ALL"
)

# B) Magnitude sensitivity: I2 rank (non-negative score).
g_i2 <- run_gsea(
  meta$ensembl_gene,
  meta$I2_total_percent,
  "pos",
  "I2_ALL"
)

# C) Coverage sensitivity: standardized-Q restricted to genes available
# in >=30 independent GSE.
m30 <- meta[N_independent_GSE >= 30]
g_q30 <- run_gsea(
  m30$ensembl_gene,
  m30$heterogeneity_Q_z,
  "std",
  "STANDARDIZED_Q_GE30_GSE"
)

fwrite(
  g_q[, !"leadingEdge", with=FALSE],
  OUT_Q, sep="\t", quote=FALSE, na="NA"
)
fwrite(
  g_i2[, !"leadingEdge", with=FALSE],
  OUT_I2, sep="\t", quote=FALSE, na="NA"
)
fwrite(
  g_q30[, !"leadingEdge", with=FALSE],
  OUT_Q30, sep="\t", quote=FALSE, na="NA"
)

# Integrate mean-response, duration-modifier, and heterogeneity pathways.
# Normalize column naming from previous analyses.
if (!"pathway" %in% names(mean_gsea))
  stop("Mean-response GSEA missing pathway")
if (!"NES" %in% names(mean_gsea))
  stop("Mean-response GSEA missing NES")
mean_fdr_col <- if ("padj" %in% names(mean_gsea)) "padj" else
  if ("padj_BH" %in% names(mean_gsea)) "padj_BH" else NA_character_
if (is.na(mean_fdr_col)) stop("Mean-response GSEA missing adjusted p-value")

if (!all(c("pathway","NES") %in% names(dur_gsea)))
  stop("Duration GSEA missing pathway/NES")
dur_fdr_col <- if ("padj_BH" %in% names(dur_gsea)) "padj_BH" else
  if ("padj" %in% names(dur_gsea)) "padj" else NA_character_
if (is.na(dur_fdr_col)) stop("Duration GSEA missing adjusted p-value")

map <- merge(
  mean_gsea[, .(
    pathway,
    mean_response_NES=NES,
    mean_response_FDR=get(mean_fdr_col)
  )],
  dur_gsea[, .(
    pathway,
    duration_modifier_NES=NES,
    duration_modifier_FDR=get(dur_fdr_col)
  )],
  by="pathway", all=TRUE
)

map <- merge(
  map,
  g_q[, .(
    pathway,
    heterogeneity_Q_NES=NES,
    heterogeneity_Q_FDR=padj_BH
  )],
  by="pathway", all=TRUE
)

map <- merge(
  map,
  g_i2[, .(
    pathway,
    heterogeneity_I2_NES=NES,
    heterogeneity_I2_FDR=padj_BH
  )],
  by="pathway", all=TRUE
)

map <- merge(
  map,
  g_q30[, .(
    pathway,
    heterogeneity_Q_ge30GSE_NES=NES,
    heterogeneity_Q_ge30GSE_FDR=padj_BH
  )],
  by="pathway", all=TRUE
)

map[, Mean_response_FDR05 := mean_response_FDR < .05]
map[, Duration_modifier_FDR05 := duration_modifier_FDR < .05]
map[, Heterogeneity_Q_FDR05 := heterogeneity_Q_FDR < .05]
map[, Heterogeneity_I2_FDR05 := heterogeneity_I2_FDR < .05]
map[, Heterogeneity_Q_ge30GSE_FDR05 := heterogeneity_Q_ge30GSE_FDR < .05]

# Descriptive labels, not inferential categories.
map[, Program_profile := fifelse(
  Mean_response_FDR05 & Heterogeneity_Q_FDR05 & Duration_modifier_FDR05,
  "MEAN_RESPONSE_PLUS_HETEROGENEITY_PLUS_DURATION",
  fifelse(
    Heterogeneity_Q_FDR05 & Duration_modifier_FDR05,
    "HETEROGENEOUS_AND_DURATION_SENSITIVE",
    fifelse(
      Mean_response_FDR05 & !Heterogeneity_Q_FDR05,
      "MEAN_RESPONSE_WITHOUT_Q_HETEROGENEITY_ENRICHMENT",
      fifelse(
        Heterogeneity_Q_FDR05,
        "HETEROGENEITY_ENRICHED",
        "OTHER"
      )
    )
  )
)]

map[, abs_het_sort := abs(heterogeneity_Q_NES)]
setorder(map, heterogeneity_Q_FDR, -abs_het_sort, pathway)
map[, abs_het_sort := NULL]
fwrite(map, OUT_MAP, sep="\t", quote=FALSE, na="NA")

fmt_top <- function(x, fdr_col="padj_BH", n=10L, positive=TRUE) {
  fdr <- x[[fdr_col]]
  if (positive) {
    z <- x[NES > 0][order(get(fdr_col), -NES)]
  } else {
    z <- x[NES < 0][order(get(fdr_col), NES)]
  }
  if (!nrow(z)) return("NONE")
  z <- head(z, n)
  paste(
    paste0(
      z$pathway,
      "|NES=", format(round(z$NES,3), nsmall=3),
      "|FDR=", format(z[[fdr_col]], digits=3, scientific=TRUE)
    ),
    collapse="\n"
  )
}

sig_q <- g_q[padj_BH < .05]
sig_i2 <- g_i2[padj_BH < .05]
sig_q30 <- g_q30[padj_BH < .05]

# Robust heterogeneity pathway = significant in standardized-Q all,
# same positive direction in I2 and >=30-GSE sensitivity, and FDR<.05
# in at least one of the two sensitivities. This is descriptive only.
rob <- merge(
  g_q[, .(pathway, NES_q=NES, FDR_q=padj_BH)],
  g_i2[, .(pathway, NES_i2=NES, FDR_i2=padj_BH)],
  by="pathway"
)
rob <- merge(
  rob,
  g_q30[, .(pathway, NES_q30=NES, FDR_q30=padj_BH)],
  by="pathway"
)
rob[, robust_heterogeneity := (
  FDR_q < .05 &
  NES_q > 0 &
  NES_i2 > 0 &
  NES_q30 > 0 &
  (FDR_i2 < .05 | FDR_q30 < .05)
)]
n_rob <- sum(rob$robust_heterogeneity, na.rm=TRUE)

summary <- c(
  "===== HFD HETEROGENEITY HALLMARK GSEA SUMMARY v1 =====",
  "",
  "STATUS=PASS",
  paste0("R_version=", R.version.string),
  paste0("fgsea_version=", as.character(packageVersion("fgsea"))),
  "MSigDB_release=2025.1.Mm",
  "MSigDB_collection=MH_Hallmark",
  "",
  "PRIMARY_EXPLORATORY_HETEROGENEITY_RANK=standardized_excess_Q",
  "STANDARDIZED_Q_FORMULA=(QE-QE_df)/sqrt(2*QE_df)",
  "INTERPRETATION_POSITIVE_NES=enrichment_toward_more_heterogeneous_genes",
  "INTERPRETATION_NEGATIVE_NES=enrichment_toward_more_homogeneous_genes",
  paste0("Ranked_genes_all=", nrow(meta)),
  paste0("Hallmark_tested_standardized_Q=", nrow(g_q)),
  paste0("Standardized_Q_FDR_lt_0.05=", nrow(sig_q)),
  paste0("Standardized_Q_positive_NES_FDR05=",
         sum(sig_q$NES > 0, na.rm=TRUE)),
  paste0("Standardized_Q_negative_NES_FDR05=",
         sum(sig_q$NES < 0, na.rm=TRUE)),
  "",
  "SENSITIVITY_1=I2_total_percent_positive_rank",
  paste0("I2_FDR_lt_0.05=", nrow(sig_i2)),
  "",
  "SENSITIVITY_2=standardized_Q_genes_ge30_independent_GSE",
  paste0("Genes_ge30_GSE=", nrow(m30)),
  paste0("Q_ge30GSE_FDR_lt_0.05=", nrow(sig_q30)),
  "",
  paste0("Descriptively_robust_high_heterogeneity_pathways=", n_rob),
  "",
  "TOP_HIGH_HETEROGENEITY_PATHWAYS_STANDARDIZED_Q",
  fmt_top(g_q, "padj_BH", 12L, TRUE),
  "",
  "TOP_LOW_HETEROGENEITY_PATHWAYS_STANDARDIZED_Q",
  fmt_top(g_q, "padj_BH", 12L, FALSE),
  "",
  paste0("Heterogeneity_rank=", sub(paste0("^",ROOT,"/"),"",OUT_RANK_Q)),
  paste0("Standardized_Q_GSEA=", sub(paste0("^",ROOT,"/"),"",OUT_Q)),
  paste0("I2_GSEA=", sub(paste0("^",ROOT,"/"),"",OUT_I2)),
  paste0("Q_ge30GSE_GSEA=", sub(paste0("^",ROOT,"/"),"",OUT_Q30)),
  paste0("Integrated_pathway_map=", sub(paste0("^",ROOT,"/"),"",OUT_MAP)),
  "",
  "ANALYSIS_STATUS=EXPLORATORY",
  "CAUSAL_INTERPRETATION=NO",
  "HUMAN_CONCORDANCE_PERFORMED=NO"
)

writeLines(summary, OUT_SUMMARY)
cat("\n", paste(summary, collapse="\n"), "\n", sep="")

