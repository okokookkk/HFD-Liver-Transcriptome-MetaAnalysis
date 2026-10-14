#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Required package not installed: data.table")
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Required package not installed: fgsea")
})

library(data.table)
library(fgsea)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

RANK_FILE <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/pathway_analysis_v1_2",
  "HFD_primary_meta_ranked_genes_zstat_v1_2.tsv.gz"
)

HALLMARK_FILE <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/pathway_analysis_v1_2",
  "MSigDB_mouse_Hallmark_frozen_v1_2.tsv.gz"
)

COVERAGE_FILE <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/pathway_analysis_v1_2",
  "MSigDB_mouse_Hallmark_rank_coverage_v1_2.tsv"
)

OUTDIR <- file.path(ROOT, "06_pathway", "hallmark_mean_response_v1_2")
METADIR <- file.path(
  ROOT,
  "00_metadata/analysis_outputs/pathway_analysis_v1"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

OUT_ALL <- file.path(
  OUTDIR,
  "HFD_Hallmark_GSEA_pooled_meta_zstat_all_v1_2.tsv"
)
OUT_SIG <- file.path(
  OUTDIR,
  "HFD_Hallmark_GSEA_pooled_meta_zstat_FDR05_v1_2.tsv"
)
OUT_LE <- file.path(
  OUTDIR,
  "HFD_Hallmark_GSEA_leading_edge_long_v1_2.tsv.gz"
)
OUT_SUMMARY <- file.path(
  METADIR,
  "HFD_Hallmark_GSEA_summary_v1_2.txt"
)

for (q in c(RANK_FILE, HALLMARK_FILE, COVERAGE_FILE)) {
  if (!file.exists(q)) stop("Missing required frozen input: ", q)
}

rank <- fread(RANK_FILE)
hall <- fread(HALLMARK_FILE)
coverage <- fread(COVERAGE_FILE)

need_rank <- c("ensembl_gene", "pathway_rank_stat")
need_hall <- c("gs_name", "ensembl_gene")
if (length(setdiff(need_rank, names(rank))))
  stop("Rank file missing required columns")
if (length(setdiff(need_hall, names(hall))))
  stop("Hallmark file missing required columns")

rank[, pathway_rank_stat := as.numeric(pathway_rank_stat)]

if (any(!is.finite(rank$pathway_rank_stat)))
  stop("Non-finite pathway rank statistic")
if (anyDuplicated(rank$ensembl_gene))
  stop("Duplicate Ensembl genes in rank file")

stats <- rank$pathway_rank_stat
names(stats) <- rank$ensembl_gene
stats <- sort(stats, decreasing = TRUE)

eligible_paths <- coverage[
  Eligible_minSize15_maxSize500 == "YES",
  gs_name
]

if (length(eligible_paths) != 50L)
  warning("Expected 50 eligible Hallmark pathways; observed ", length(eligible_paths))

hall_use <- hall[gs_name %in% eligible_paths]
pathways <- split(hall_use$ensembl_gene, hall_use$gs_name)

# Deterministic seed for auditability.
set.seed(20260904)

cores_env <- Sys.getenv("HFD_FGSEA_CORES", unset = "")
if (nzchar(cores_env)) {
  NCORES <- as.integer(cores_env)
} else {
  NCORES <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
}
if (!is.finite(NCORES) || NCORES < 1L) NCORES <- 1L

cat("===== HFD HALLMARK GSEA v1.2 =====\n")
cat("fgsea:", as.character(packageVersion("fgsea")), "\n")
cat("Pathways:", length(pathways), "\n")
cat("Ranked genes:", length(stats), "\n")
cat("Workers:", NCORES, "\n\n")

fg <- fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  scoreType = "std",
  nproc = NCORES,
  gseaParam = 1
)

res <- as.data.table(fg)

if (nrow(res) != length(pathways))
  warning("fgsea returned ", nrow(res), " pathways for ", length(pathways), " inputs")

# Attach frozen overlap count.
res <- merge(
  res,
  coverage[, .(
    pathway = gs_name,
    N_MSigDB_genes,
    N_ranked_genes_overlap,
    Overlap_fraction
  )],
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)

# Explicit BH correction across the 50 Hallmark pathways.
res[, padj_BH_50 := p.adjust(pval, method = "BH")]

# fgsea's padj should match BH aside from any internal NA handling.
padj_diff <- max(abs(res$padj - res$padj_BH_50), na.rm = TRUE)

# Collapse leading edge for the compact all-pathway table.
res[, leadingEdge_ensembl := vapply(
  leadingEdge,
  function(z) paste(z, collapse = ";"),
  character(1)
)]
res[, leadingEdge_size := lengths(leadingEdge)]

# Keep fgsea output plus explicit audit fields.
out <- res[, .(
  pathway,
  NES,
  ES,
  pval,
  padj = padj_BH_50,
  log2err,
  size,
  N_MSigDB_genes,
  N_ranked_genes_overlap,
  Overlap_fraction,
  leadingEdge_size,
  leadingEdge_ensembl
)]

out[, abs_NES_sort := abs(NES)]
setorder(out, padj, -abs_NES_sort, pathway)
out[, abs_NES_sort := NULL]
fwrite(out, OUT_ALL, sep = "\t", quote = FALSE, na = "NA")

sig <- out[padj < 0.05]
sig[, abs_NES_sort := abs(NES)]
setorder(sig, padj, -abs_NES_sort, pathway)
sig[, abs_NES_sort := NULL]
fwrite(sig, OUT_SIG, sep = "\t", quote = FALSE, na = "NA")

# Long-form leading edge with rank statistics and annotations.
if (nrow(sig)) {
  le <- res[padj_BH_50 < 0.05, .(
    ensembl_gene = unlist(leadingEdge)
  ), by = .(pathway, NES, pval, padj = padj_BH_50)]

  le <- merge(
    le,
    rank[, .(
      ensembl_gene,
      Geneid,
      gene_symbol,
      gene_type,
      pathway_rank_stat,
      pooled_log2FC,
      pooled_SE,
      pooled_FDR_BH,
      N_independent_GSE,
      I2_total_percent
    )],
    by = "ensembl_gene",
    all.x = TRUE,
    sort = FALSE
  )

  le[, abs_pathway_rank_stat_sort := abs(pathway_rank_stat)]
  setorder(le, pathway, -abs_pathway_rank_stat_sort, ensembl_gene)
  le[, abs_pathway_rank_stat_sort := NULL]
  fwrite(le, OUT_LE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
} else {
  empty <- data.table(
    pathway = character(),
    NES = numeric(),
    pval = numeric(),
    padj = numeric(),
    ensembl_gene = character(),
    Geneid = character(),
    gene_symbol = character(),
    gene_type = character(),
    pathway_rank_stat = numeric(),
    pooled_log2FC = numeric(),
    pooled_SE = numeric(),
    pooled_FDR_BH = numeric(),
    N_independent_GSE = integer(),
    I2_total_percent = numeric()
  )
  fwrite(empty, OUT_LE, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
}

n_sig <- nrow(sig)
n_pos <- sum(sig$NES > 0, na.rm = TRUE)
n_neg <- sum(sig$NES < 0, na.rm = TRUE)
n_fdr01 <- sum(out$padj < 0.01, na.rm = TRUE)
n_nom <- sum(out$pval < 0.05, na.rm = TRUE)

top_pos <- out[NES > 0][order(padj, -NES)]
top_neg <- out[NES < 0][order(padj, NES)]

fmt_top <- function(x, n = 8L) {
  if (nrow(x) == 0L) return("NONE")
  x <- head(x, n)
  paste(
    paste0(
      x$pathway,
      "|NES=", format(round(x$NES, 3), nsmall = 3),
      "|FDR=", format(x$padj, digits = 3, scientific = TRUE)
    ),
    collapse = "\n"
  )
}

summary <- c(
  "===== HFD HALLMARK GSEA SUMMARY v1.2 =====",
  "",
  "STATUS=PASS",
  paste0("R_version=", R.version.string),
  paste0("fgsea_version=", as.character(packageVersion("fgsea"))),
  paste0("Random_seed=20260904"),
  paste0("Parallel_workers=", NCORES),
  "",
  "RANKING=pooled_meta_z_statistic",
  "ENGINE=fgseaMultilevel",
  "SCORE_TYPE=std",
  "GSEA_PARAM=1",
  "MIN_SIZE=15",
  "MAX_SIZE=500",
  "",
  paste0("Hallmark_pathways_tested=", nrow(out)),
  paste0("Nominal_p_lt_0.05=", n_nom),
  paste0("BH_FDR_lt_0.05=", n_sig),
  paste0("BH_FDR_lt_0.01=", n_fdr01),
  paste0("FDR05_positive_NES=", n_pos),
  paste0("FDR05_negative_NES=", n_neg),
  paste0("Max_abs_fgsea_vs_explicit_BH_padj_difference=",
         format(padj_diff, digits = 8)),
  "",
  "TOP_POSITIVE_BY_FDR",
  fmt_top(top_pos, 8L),
  "",
  "TOP_NEGATIVE_BY_FDR",
  fmt_top(top_neg, 8L),
  "",
  paste0("All_results=", sub(paste0("^", ROOT, "/"), "", OUT_ALL)),
  paste0("FDR05_results=", sub(paste0("^", ROOT, "/"), "", OUT_SIG)),
  paste0("Leading_edge_long=", sub(paste0("^", ROOT, "/"), "", OUT_LE)),
  "",
  "HETEROGENEITY_PATHWAY_ANALYSIS_PERFORMED=NO",
  "MODERATOR_PATHWAY_ANALYSIS_PERFORMED=NO"
)

writeLines(summary, OUT_SUMMARY)
cat(paste(summary, collapse = "\n"), "\n")

