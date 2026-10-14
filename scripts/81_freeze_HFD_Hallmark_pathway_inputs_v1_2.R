#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Required R package not installed: data.table")
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Required R package not installed: msigdbr")
})

library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

META <- file.path(
  ROOT, "05_meta/gene_level_v1_2",
  "HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz"
)

OUTDIR <- file.path(
  ROOT, "00_metadata/analysis_inputs/pathway_analysis_v1_2"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUT_RANK <- file.path(
  OUTDIR, "HFD_primary_meta_ranked_genes_zstat_v1_2.tsv.gz"
)
OUT_HALLMARK <- file.path(
  OUTDIR, "MSigDB_mouse_Hallmark_frozen_v1_2.tsv.gz"
)
OUT_COVERAGE <- file.path(
  OUTDIR, "MSigDB_mouse_Hallmark_rank_coverage_v1_2.tsv"
)
OUT_SUMMARY <- file.path(
  OUTDIR, "HFD_Hallmark_pathway_input_freeze_summary_v1_2.txt"
)

if (!file.exists(META)) stop("Missing primary meta-analysis result: ", META)

meta <- fread(META)

need <- c(
  "Geneid", "gene_symbol", "gene_type",
  "pooled_log2FC", "pooled_SE", "z_value", "p_value", "pooled_FDR_BH",
  "N_independent_GSE", "tau2_total", "I2_total_percent",
  "prediction95_lower_new_GSE", "prediction95_upper_new_GSE"
)
missing <- setdiff(need, names(meta))
if (length(missing))
  stop("Meta table missing columns: ", paste(missing, collapse = ", "))

if (nrow(meta) != 18234L)
  warning("Expected 18,234 primary meta genes; observed ", nrow(meta))

# ------------------------------------------------------------------
# Freeze genome-wide signed ranking.
# Primary ranking statistic: pooled meta-analysis z statistic.
# ------------------------------------------------------------------
rank <- meta[, ..need]
rank[, ensembl_gene := sub("\\.[0-9]+$", "", Geneid)]

if (anyDuplicated(rank$ensembl_gene))
  stop("Version-stripped Ensembl gene IDs are not unique in primary meta table")

numcols <- c(
  "pooled_log2FC", "pooled_SE", "z_value", "p_value", "pooled_FDR_BH",
  "N_independent_GSE", "tau2_total", "I2_total_percent",
  "prediction95_lower_new_GSE", "prediction95_upper_new_GSE"
)
for (cc in numcols) set(rank, j = cc, value = as.numeric(rank[[cc]]))

if (any(!is.finite(rank$z_value)))
  stop("Non-finite z statistic in primary meta ranking")
if (any(!is.finite(rank$pooled_log2FC)) || any(!is.finite(rank$pooled_SE)))
  stop("Non-finite pooled effect/SE in primary meta table")
if (any(rank$pooled_SE <= 0))
  stop("Non-positive pooled SE")

z_check <- rank$pooled_log2FC / rank$pooled_SE
max_z_diff <- max(abs(rank$z_value - z_check))
if (!is.finite(max_z_diff) || max_z_diff > 1e-6)
  stop("Stored z_value differs from pooled_log2FC/pooled_SE; max diff=", max_z_diff)

rank[, pathway_rank_stat := z_value]
setorder(rank, -pathway_rank_stat, ensembl_gene)

rank <- rank[, .(
  Geneid,
  ensembl_gene,
  gene_symbol,
  gene_type,
  pathway_rank_stat,
  pooled_log2FC,
  pooled_SE,
  z_value,
  p_value,
  pooled_FDR_BH,
  N_independent_GSE,
  tau2_total,
  I2_total_percent,
  prediction95_lower_new_GSE,
  prediction95_upper_new_GSE
)]

fwrite(rank, OUT_RANK, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

# ------------------------------------------------------------------
# Freeze native mouse MSigDB Hallmark collection from local msigdbr.
# IMPORTANT: native mouse collection codes differ from human codes
# (e.g. mouse Hallmark is currently MH, human Hallmark is H).
# Therefore discover the Hallmark collection by collection NAME rather
# than hard-coding a collection abbreviation.
# Do NOT silently substitute human-ortholog Hallmark.
# ------------------------------------------------------------------
collections <- tryCatch(
  msigdbr::msigdbr_collections(db_species = "MM"),
  error = function(e) e
)

if (inherits(collections, "error")) {
  stop(
    "Could not query native mouse MSigDB collections via msigdbr_collections(). ",
    "Original error: ", conditionMessage(collections)
  )
}

collections <- as.data.table(collections)

# Older msigdbr releases may not expose db_version from
# msigdbr_collections(). Collection discovery only requires these two
# fields; the actual database version is frozen later from msigdbr().
required_collections <- c(
  "gs_collection", "gs_collection_name"
)
missing_collections <- setdiff(required_collections, names(collections))
if (length(missing_collections))
  stop(
    "msigdbr_collections output missing required collection columns: ",
    paste(missing_collections, collapse = ", "),
    "; observed columns: ", paste(names(collections), collapse = ", ")
  )

hallmark_collection <- unique(
  collections[gs_collection_name == "Hallmark", gs_collection]
)

if (length(hallmark_collection) != 1L) {
  stop(
    "Expected exactly one native mouse Hallmark collection; observed: ",
    paste(hallmark_collection, collapse = ", ")
  )
}

hallmark_collection_name <- "Hallmark"

msig <- tryCatch(
  msigdbr::msigdbr(
    db_species = "MM",
    species = "Mus musculus",
    collection = hallmark_collection
  ),
  error = function(e) e
)

if (inherits(msig, "error")) {
  stop(
    "Installed msigdbr could not retrieve native mouse Hallmark collection ",
    hallmark_collection, ". Do not substitute human orthologs automatically. ",
    "Original error: ", conditionMessage(msig)
  )
}

msig <- as.data.table(msig)

required_msig <- c(
  "gene_symbol", "ensembl_gene", "gs_name",
  "gs_collection"
)
missing_msig <- setdiff(required_msig, names(msig))
if (length(missing_msig))
  stop(
    "msigdbr output missing required columns: ",
    paste(missing_msig, collapse = ", "),
    "; observed columns: ", paste(names(msig), collapse = ", ")
  )

if (!"db_version" %in% names(msig)) {
  stop(
    "Installed msigdbr retrieved Hallmark sets but does not expose db_version. ",
    "Exact MSigDB release cannot be frozen safely. ",
    "msigdbr package version=", as.character(packageVersion("msigdbr")),
    "; observed columns: ", paste(names(msig), collapse = ", ")
  )
}

hall <- msig[gs_collection == hallmark_collection]

if (nrow(hall) == 0L)
  stop(
    "No Hallmark records returned from native mouse MSigDB collection ",
    hallmark_collection
  )

versions <- unique(hall$db_version)
versions <- versions[!is.na(versions) & nzchar(versions)]
if (length(versions) != 1L)
  stop("Expected exactly one MSigDB db_version for Hallmark; observed: ",
       paste(versions, collapse = ", "))

keep_optional <- intersect(
  c("gs_description", "gs_source_species", "db_target_species"),
  names(hall)
)

hall <- hall[, c(
  "db_version", "gs_name", keep_optional,
  "gene_symbol", "ensembl_gene"
), with = FALSE]

hall <- hall[
  !is.na(ensembl_gene) & nzchar(ensembl_gene) &
  !is.na(gs_name) & nzchar(gs_name)
]

hall <- unique(hall, by = c("gs_name", "ensembl_gene"))
setorder(hall, gs_name, ensembl_gene)

fwrite(hall, OUT_HALLMARK, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

# ------------------------------------------------------------------
# Audit overlap with the frozen ranked gene universe.
# ------------------------------------------------------------------
rank_ids <- unique(rank$ensembl_gene)

coverage <- hall[, .(
  N_MSigDB_genes = uniqueN(ensembl_gene),
  N_ranked_genes_overlap = uniqueN(ensembl_gene[ensembl_gene %in% rank_ids])
), by = gs_name]

coverage[, Overlap_fraction := N_ranked_genes_overlap / N_MSigDB_genes]
coverage[, Eligible_minSize15_maxSize500 := ifelse(
  N_ranked_genes_overlap >= 15 & N_ranked_genes_overlap <= 500,
  "YES", "NO"
)]
setorder(coverage, gs_name)
fwrite(coverage, OUT_COVERAGE, sep = "\t", quote = FALSE, na = "NA")

n_pathways <- uniqueN(hall$gs_name)
n_eligible <- sum(coverage$Eligible_minSize15_maxSize500 == "YES")
n_rank_pc <- sum(rank$gene_type == "protein_coding", na.rm = TRUE)
n_rank_nonpc <- nrow(rank) - n_rank_pc
n_rank_msig <- sum(rank$ensembl_gene %in% unique(hall$ensembl_gene))

summary <- c(
  "===== HFD HALLMARK PATHWAY INPUT FREEZE SUMMARY v1.2 =====",
  "",
  "STATUS=PASS",
  paste0("R_version=", R.version.string),
  paste0("msigdbr_version=", as.character(packageVersion("msigdbr"))),
  paste0("MSigDB_database_species=MM"),
  paste0("MSigDB_target_species=Mus musculus"),
  paste0("MSigDB_collection_code=", hallmark_collection),
  paste0("MSigDB_collection_name=", hallmark_collection_name),
  paste0("MSigDB_db_version=", versions),
  "",
  "PRIMARY_PATHWAY_RANKING=pooled_meta_z_statistic",
  "RANK_FORMULA=pooled_log2FC_divided_by_pooled_SE",
  paste0("Ranked_genes=", nrow(rank)),
  paste0("Ranked_protein_coding_genes=", n_rank_pc),
  paste0("Ranked_non_protein_coding_genes=", n_rank_nonpc),
  paste0("Ranked_genes_overlapping_Hallmark=", n_rank_msig),
  paste0("Max_abs_z_recalculation_difference=", format(max_z_diff, digits = 8)),
  "",
  paste0("Hallmark_pathways_total=", n_pathways),
  paste0("Hallmark_pathways_eligible_min15_max500=", n_eligible),
  paste0("Minimum_Hallmark_overlap=", min(coverage$N_ranked_genes_overlap)),
  paste0("Median_Hallmark_overlap=", median(coverage$N_ranked_genes_overlap)),
  paste0("Maximum_Hallmark_overlap=", max(coverage$N_ranked_genes_overlap)),
  "",
  "PRIMARY_GSEA_MIN_SIZE=15",
  "PRIMARY_GSEA_MAX_SIZE=500",
  "PRIMARY_GSEA_SCORE_TYPE=std",
  "PRIMARY_GSEA_GSEA_PARAM=1",
  "PRIMARY_GSEA_MULTIPLE_TESTING=BH_across_Hallmark_pathways",
  "PRIMARY_GSEA_ENGINE=fgseaMultilevel",
  "",
  "HUMAN_ORTHOLOG_GENE_SETS_USED=NO",
  "STABLE_CORE_59_USED_AS_PRIMARY_PATHWAY_INPUT=NO",
  "PATHWAY_ENRICHMENT_PERFORMED=NO",
  "",
  paste0("Frozen_ranked_genes=", sub(paste0("^", ROOT, "/"), "", OUT_RANK)),
  paste0("Frozen_Hallmark_gene_sets=", sub(paste0("^", ROOT, "/"), "", OUT_HALLMARK)),
  paste0("Hallmark_rank_coverage=", sub(paste0("^", ROOT, "/"), "", OUT_COVERAGE))
)

writeLines(summary, OUT_SUMMARY)
cat(paste(summary, collapse = "\n"), "\n")

