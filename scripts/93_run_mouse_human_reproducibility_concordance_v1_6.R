#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table","AnnotationDbi","org.Hs.eg.db")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse=", "))
})
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

MOUSE <- file.path(
  ROOT, "05_meta/gene_level_v1_2",
  "HFD_gene_multilevel_meta_primary_annotated_v1_2.tsv.gz"
)
MOUSE_CORE <- file.path(
  ROOT, "05_meta/stable_core_v1",
  "HFD_stable_core_final_v1.tsv"
)
HUMAN <- file.path(
  ROOT, "07_human_validation/gene_meta_v1",
  "human_MASLD_three_cohort_random_effects_meta_v1.tsv.gz"
)
HUMAN_STRICT <- file.path(
  ROOT, "07_human_validation/gene_meta_v1",
  "human_MASLD_three_cohort_strict4_sensitivity_meta_v1.tsv.gz"
)

OUTDIR <- file.path(ROOT, "07_human_validation", "mouse_human_concordance_v1_6")
METADIR <- file.path(
  ROOT, "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_v1_6"
)
MAPDIR <- file.path(
  ROOT, "00_metadata/human_validation/analysis_inputs/ortholog_map_v1_6"
)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)
dir.create(MAPDIR, recursive=TRUE, showWarnings=FALSE)

ORTHO <- file.path(
  MAPDIR, "Ensembl116_mouse_human_highconfidence_one2one_orthologs_v1_6.tsv.gz"
)
JOINED <- file.path(
  OUTDIR, "mouse_HFD_human_MASLD_one2one_ortholog_concordance_v1_6.tsv.gz"
)
QUART <- file.path(
  OUTDIR, "mouse_reproducibility_quartile_human_concordance_v1_6.tsv"
)
MODELS <- file.path(
  OUTDIR, "mouse_reproducibility_predicts_human_concordance_models_v1_6.tsv"
)
COREOUT <- file.path(
  OUTDIR, "mouse_stable_core_human_concordance_descriptive_v1_6.tsv"
)
SUMMARY <- file.path(
  METADIR, "mouse_human_reproducibility_concordance_summary_v1_6.txt"
)

for (q in c(MOUSE,HUMAN,HUMAN_STRICT)) {
  if (!file.exists(q)) stop("Missing required input: ", q)
}

mouse <- fread(MOUSE)
human <- fread(HUMAN)
humanS <- fread(HUMAN_STRICT)

# -------------------------------------------------------------------------
# Robust column resolver. Stop rather than guess if a required field cannot
# be uniquely identified.
# -------------------------------------------------------------------------
resolve_col <- function(dt, exact, regex=NULL, label, required=TRUE) {
  hit <- exact[exact %in% names(dt)]
  if (length(hit) == 1L) return(hit)
  if (length(hit) > 1L) return(hit[1L])

  if (!is.null(regex)) {
    rhit <- grep(regex, names(dt), value=TRUE, ignore.case=TRUE, perl=TRUE)
    if (length(rhit) == 1L) return(rhit)
    if (length(rhit) > 1L) {
      stop(
        "Ambiguous column resolution for ", label, ": ",
        paste(rhit, collapse=", "),
        "\nAvailable columns:\n", paste(names(dt), collapse="\n")
      )
    }
  }

  if (required) {
    stop(
      "Could not resolve required column ", label,
      "\nAvailable columns:\n", paste(names(dt), collapse="\n")
    )
  }
  NA_character_
}

m_gene <- resolve_col(
  mouse,
  c("gene_id","Geneid","ensembl_gene_id","ENSEMBL","ENSEMBL_gene_id"),
  "^(gene_?id|geneid|ensembl.*gene.*id)$",
  "mouse Ensembl gene ID"
)
m_type <- resolve_col(
  mouse,
  c("gene_type","gene_biotype","Gene_type","biotype"),
  "(gene.*type|gene.*biotype|^biotype$)",
  "mouse gene type"
)
m_beta <- resolve_col(
  mouse,
  c("pooled_log2FC","pooled_effect","meta_log2FC","pooled_beta","estimate"),
  "(pooled.*log2fc|meta.*log2fc|pooled.*effect|pooled.*beta)",
  "mouse pooled log2FC"
)
m_se <- resolve_col(
  mouse,
  c("pooled_SE","meta_SE","pooled_se"),
  "(pooled|meta).*[_\\.](se)$|^(pooled|meta).*se$",
  "mouse pooled SE"
)
m_ngse <- resolve_col(
  mouse,
  c(
    "n_GSE","N_GSE","n_independent_GSE","independent_GSE",
    "N_independent_GSE","k_GSE","n_studies"
  ),
  "(n|num|number|k).*(independent.*)?gse|independent.*gse.*(n|count)",
  "mouse independent GSE count"
)
m_dir <- resolve_col(
  mouse,
  c(
    "GSE_balanced_direction_consistency",
    "gse_balanced_direction_consistency",
    "GSE_direction_consistency",
    "gse_direction_consistency",
    "direction_consistency_GSE"
  ),
  "(gse.*direction.*consistency|direction.*consistency.*gse)",
  "mouse GSE-balanced direction consistency"
)
m_i2 <- resolve_col(
  mouse,
  c("I2_total_percent","I2_total","i2_total_percent","i2_total"),
  "(^|_)i2.*total|total.*i2",
  "mouse total I2"
)

m_pi_lo <- resolve_col(
  mouse,
  c(
    "prediction_interval_lower","PI_lower","pred_interval_lower",
    "prediction_interval_lb","PI_lb","pred_lb","pi_lower"
  ),
  "(prediction|pred|pi).*(lower|lb)",
  "mouse PI lower", required=FALSE
)
m_pi_hi <- resolve_col(
  mouse,
  c(
    "prediction_interval_upper","PI_upper","pred_interval_upper",
    "prediction_interval_ub","PI_ub","pred_ub","pi_upper"
  ),
  "(prediction|pred|pi).*(upper|ub)",
  "mouse PI upper", required=FALSE
)

cat("===== MOUSE-HUMAN REPRODUCIBILITY CONCORDANCE v1.6 =====\n")
cat("Resolved mouse columns:\n")
cat(" gene=",m_gene,"\n type=",m_type,"\n beta=",m_beta,"\n se=",m_se,
    "\n nGSE=",m_ngse,"\n direction=",m_dir,"\n I2=",m_i2,"\n",sep="")

# -------------------------------------------------------------------------
# Mouse protein-coding primary-meta universe.
# -------------------------------------------------------------------------
m <- mouse[, .(
  mouse_ensembl=sub("\\.[0-9]+$","",as.character(get(m_gene))),
  mouse_gene_type=as.character(get(m_type)),
  mouse_log2FC=as.numeric(get(m_beta)),
  mouse_SE=as.numeric(get(m_se)),
  mouse_nGSE=as.numeric(get(m_ngse)),
  mouse_direction_consistency=as.numeric(get(m_dir)),
  mouse_I2=as.numeric(get(m_i2))
)]

if (!is.na(m_pi_lo) && !is.na(m_pi_hi)) {
  m[, mouse_PI_lower := as.numeric(mouse[[m_pi_lo]])]
  m[, mouse_PI_upper := as.numeric(mouse[[m_pi_hi]])]
} else {
  m[, mouse_PI_lower := NA_real_]
  m[, mouse_PI_upper := NA_real_]
}

m <- m[
  mouse_gene_type == "protein_coding" &
  nzchar(mouse_ensembl) &
  is.finite(mouse_log2FC) &
  is.finite(mouse_SE) & mouse_SE > 0 &
  is.finite(mouse_nGSE) &
  is.finite(mouse_direction_consistency) &
  is.finite(mouse_I2)
]

if (anyDuplicated(m$mouse_ensembl))
  stop("Duplicate protein-coding mouse Ensembl IDs after annotation filtering")

# Normalize direction consistency if stored as percent.
if (max(m$mouse_direction_consistency, na.rm=TRUE) > 1.5) {
  m[, mouse_direction_consistency := mouse_direction_consistency / 100]
}
if (any(m$mouse_direction_consistency < 0 | m$mouse_direction_consistency > 1))
  stop("Mouse direction consistency outside [0,1] after normalization")

# Normalize I2 if stored as fraction.
if (max(m$mouse_I2, na.rm=TRUE) <= 1.5) {
  m[, mouse_I2_percent := 100 * mouse_I2]
} else {
  m[, mouse_I2_percent := mouse_I2]
}

m[, mouse_z := mouse_log2FC / mouse_SE]
m[, mouse_abs_log2FC := abs(mouse_log2FC)]
m[, mouse_PI_excludes_zero := (
  is.finite(mouse_PI_lower) & is.finite(mouse_PI_upper) &
  ((mouse_PI_lower > 0 & mouse_log2FC > 0) |
   (mouse_PI_upper < 0 & mouse_log2FC < 0))
)]

# -------------------------------------------------------------------------
# Reuse the already generated Ensembl release-116 high-confidence one-to-one
# ortholog map from v1_4.
#
# v1_4 successfully downloaded and filtered the official Ensembl 116 Compara
# FTP homology dumps and wrote the frozen ortholog table. Its only failure was
# later, during summary construction, because a stale ENSEMBL_HOST variable
# remained. To avoid unnecessary network dependence, v1_5 reuses only that
# validated frozen ortholog map and reruns the concordance analysis.
# -------------------------------------------------------------------------
ENSEMBL_RELEASE <- 116L

ORTHO_V1_4 <- file.path(
  ROOT,
  "00_metadata/human_validation/analysis_inputs/ortholog_map_v1_4",
  "Ensembl116_mouse_human_highconfidence_one2one_orthologs_v1_4.tsv.gz"
)

if (!file.exists(ORTHO_V1_4))
  stop("Missing successfully generated v1_4 frozen ortholog map: ", ORTHO_V1_4)

orth <- fread(ORTHO_V1_4)

# Canonical join-key typing: IDs are identifiers, never numeric quantities.
# fread can infer ENTREZID as integer; force all cross-dataset IDs to character.
orth[, mouse_ensembl := as.character(mouse_ensembl)]
orth[, human_ensembl := as.character(human_ensembl)]
orth[, human_ENTREZID := as.character(human_ENTREZID)]

orth_required <- c(
  "mouse_ensembl","mouse_symbol",
  "human_ensembl","human_ENTREZID","human_symbol",
  "orthology_type","orthology_confidence"
)
orth_missing <- setdiff(orth_required, names(orth))
if (length(orth_missing))
  stop(
    "v1_4 frozen ortholog map missing columns: ",
    paste(orth_missing, collapse=",")
  )

orth <- orth[, ..orth_required]

if (nrow(orth) < 10000L)
  stop("Unexpectedly few v1_4 frozen one-to-one orthologs: ", nrow(orth))

if (any(
  orth$orthology_type != "ortholog_one2one" |
  as.character(orth$orthology_confidence) != "1"
))
  stop("v1_4 frozen ortholog map violates one-to-one/high-confidence rule")

if (anyDuplicated(orth$mouse_ensembl))
  stop("Duplicate mouse Ensembl IDs in v1_4 frozen ortholog map")
if (anyDuplicated(orth$human_ensembl))
  stop("Duplicate human Ensembl IDs in v1_4 frozen ortholog map")
if (anyDuplicated(orth$human_ENTREZID))
  stop("Duplicate human Entrez IDs in v1_4 frozen ortholog map")

# Freeze an explicit v1_5 copy for a self-contained audit trail.
fwrite(orth, ORTHO, sep="\t", quote=FALSE, na="NA", compress="gzip")

orth_source_sha256 <- system2(
  "sha256sum",
  shQuote(ORTHO_V1_4),
  stdout=TRUE
)
orth_v1_5_sha256 <- system2(
  "sha256sum",
  shQuote(ORTHO),
  stdout=TRUE
)
orth_source_sha256 <- sub("\\s+.*$", "", orth_source_sha256[1])
orth_v1_5_sha256 <- sub("\\s+.*$", "", orth_v1_5_sha256[1])

if (!identical(orth_source_sha256, orth_v1_5_sha256))
  stop("v1_5 ortholog copy SHA256 differs from v1_4 frozen source")

# -------------------------------------------------------------------------
# Join mouse and human primary + strict4 effects.
# Canonicalize all join keys to character before any merge.
m[, mouse_ensembl := as.character(mouse_ensembl)]
orth[, mouse_ensembl := as.character(mouse_ensembl)]
orth[, human_ENTREZID := as.character(human_ENTREZID)]

# -------------------------------------------------------------------------
h <- human[, .(
  human_ENTREZID=as.character(ENTREZID),
  human_RE_log2FC=as.numeric(human_RE_log2FC),
  human_RE_SE=as.numeric(human_RE_SE),
  human_RE_FDR=as.numeric(human_RE_FDR_BH),
  human_FE_log2FC=as.numeric(human_FE_log2FC),
  human_FE_FDR=as.numeric(human_FE_FDR_BH),
  human_all3_same=as.logical(all_three_same_direction),
  human_all_LOO_FE_same_sign=as.logical(all_LOO_FE_same_sign)
)]
hs <- humanS[, .(
  human_ENTREZID=as.character(ENTREZID),
  human_strict4_RE_log2FC=as.numeric(human_RE_log2FC),
  human_strict4_RE_FDR=as.numeric(human_RE_FDR_BH)
)]

if (!is.character(m$mouse_ensembl))
  stop("mouse_ensembl is not character before join")
if (!is.character(orth$mouse_ensembl))
  stop("orth mouse_ensembl is not character before join")
if (!is.character(orth$human_ENTREZID))
  stop("orth human_ENTREZID is not character before join")
if (!is.character(h$human_ENTREZID))
  stop("human primary ENTREZID is not character before join")
if (!is.character(hs$human_ENTREZID))
  stop("human strict4 ENTREZID is not character before join")

z <- merge(m,orth,by="mouse_ensembl")
z <- merge(z,h,by="human_ENTREZID")
z <- merge(z,hs,by="human_ENTREZID",all.x=TRUE)

if (!nrow(z)) stop("No mouse-human orthologs survived the joined analysis universe")
if (anyDuplicated(z$mouse_ensembl) || anyDuplicated(z$human_ENTREZID))
  stop("Duplicate orthologs after mouse-human join")

z[, concordant_primary := sign(mouse_log2FC) == sign(human_RE_log2FC)]
z[, concordant_FE := sign(mouse_log2FC) == sign(human_FE_log2FC)]
z[, concordant_strict4 := sign(mouse_log2FC) == sign(human_strict4_RE_log2FC)]
z[, human_aligned_effect := sign(mouse_log2FC) * human_RE_log2FC]
z[, human_aligned_FE_effect := sign(mouse_log2FC) * human_FE_log2FC]
z[, human_aligned_strict4_effect := sign(mouse_log2FC) * human_strict4_RE_log2FC]
z[, human_significant_concordant := concordant_primary & human_RE_FDR < 0.05]
z[, human_all3_concordant := human_all3_same & concordant_primary]

# Prespecified primary predictor scale: per +10 percentage points direction consistency.
z[, mouse_direction_consistency_10pp := 10 * mouse_direction_consistency]
z[, mouse_I2_10pp := mouse_I2_percent / 10]
z[, log10_mouse_nGSE := log10(mouse_nGSE)]

# -------------------------------------------------------------------------
# Stable-core flag is secondary/descriptive only.
# -------------------------------------------------------------------------
z[, mouse_stable_core := FALSE]
if (file.exists(MOUSE_CORE)) {
  core <- fread(MOUSE_CORE)
  c_gene <- resolve_col(
    core,
    c("gene_id","Geneid","ensembl_gene_id","ENSEMBL","ENSEMBL_gene_id"),
    "^(gene_?id|geneid|ensembl.*gene.*id)$",
    "stable-core mouse Ensembl gene ID"
  )
  core_ids <- unique(sub("\\.[0-9]+$","",as.character(core[[c_gene]])))
  z[, mouse_stable_core := mouse_ensembl %in% core_ids]
}

# -------------------------------------------------------------------------
# Primary model:
# Human RE same-direction ~ mouse GSE-balanced direction consistency.
# Adjusted model asks whether this persists beyond mouse effect magnitude
# and study coverage. I2 is secondary, not entered with direction consistency.
# -------------------------------------------------------------------------
fit_glm <- function(formula, dt, label, predictor) {
  fit <- glm(formula, data=dt, family=binomial())
  sm <- summary(fit)$coefficients
  if (!(predictor %in% rownames(sm))) stop(label, ": predictor absent in model")
  data.table(
    Model=label,
    Predictor=predictor,
    N=nobs(fit),
    Beta_logOR=sm[predictor,"Estimate"],
    SE=sm[predictor,"Std. Error"],
    Z=sm[predictor,"z value"],
    P=sm[predictor,"Pr(>|z|)"],
    OR=exp(sm[predictor,"Estimate"]),
    OR_95CI_low=exp(sm[predictor,"Estimate"] - 1.96*sm[predictor,"Std. Error"]),
    OR_95CI_high=exp(sm[predictor,"Estimate"] + 1.96*sm[predictor,"Std. Error"])
  )
}

models <- rbindlist(list(
  fit_glm(
    concordant_primary ~ mouse_direction_consistency_10pp,
    z,
    "PRIMARY_unadjusted_RE_same_direction",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    concordant_primary ~ mouse_direction_consistency_10pp +
      mouse_abs_log2FC + log10_mouse_nGSE,
    z,
    "PRIMARY_adjusted_for_mouse_effect_and_coverage",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    concordant_primary ~ mouse_I2_10pp,
    z,
    "SECONDARY_mouse_I2",
    "mouse_I2_10pp"
  ),
  fit_glm(
    concordant_FE ~ mouse_direction_consistency_10pp,
    z,
    "SENSITIVITY_human_FE_direction",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    concordant_strict4 ~ mouse_direction_consistency_10pp,
    z[is.finite(human_strict4_RE_log2FC)],
    "SENSITIVITY_GSE130970_strict4_human_RE_direction",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    concordant_primary ~ mouse_direction_consistency_10pp,
    z[mouse_nGSE >= 20],
    "SENSITIVITY_mouse_coverage_ge20_GSE",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    human_significant_concordant ~ mouse_direction_consistency_10pp,
    z,
    "SECONDARY_human_FDR05_and_same_direction",
    "mouse_direction_consistency_10pp"
  ),
  fit_glm(
    human_all3_concordant ~ mouse_direction_consistency_10pp,
    z,
    "SECONDARY_all3_human_cohorts_same_mouse_direction",
    "mouse_direction_consistency_10pp"
  )
), fill=TRUE)
models[, FDR_BH_secondary_family := p.adjust(P, method="BH")]
fwrite(models, MODELS, sep="\t", quote=FALSE, na="NA")

# -------------------------------------------------------------------------
# Reproducibility quartiles for transparent effect-size presentation.
# Quartiles are based on the joined ortholog universe.
# -------------------------------------------------------------------------
# Rank-balanced quartiles avoid duplicated quantile-break errors when
# direction-consistency values are highly discrete/tied.
z[, repro_rank := frank(
  mouse_direction_consistency,
  ties.method="average",
  na.last="keep"
)]
z[, repro_q_index := pmin(
  4L,
  pmax(1L, ceiling(4 * repro_rank / .N))
)]
z[, reproducibility_quartile := factor(
  c("Q1_low","Q2","Q3","Q4_high")[repro_q_index],
  levels=c("Q1_low","Q2","Q3","Q4_high")
)]
z[, c("repro_rank","repro_q_index") := NULL]

quart <- z[, .(
  N=.N,
  Median_mouse_direction_consistency=median(mouse_direction_consistency),
  Median_mouse_nGSE=median(mouse_nGSE),
  Human_same_direction_fraction=mean(concordant_primary),
  Human_FDR05_same_direction_fraction=mean(human_significant_concordant),
  Human_all3_concordant_fraction=mean(human_all3_concordant),
  Median_human_aligned_effect=median(human_aligned_effect),
  Mean_human_aligned_effect=mean(human_aligned_effect)
), by=reproducibility_quartile]
fwrite(quart, QUART, sep="\t", quote=FALSE, na="NA")

# -------------------------------------------------------------------------
# Continuous descriptive associations.
# -------------------------------------------------------------------------
sp_dir_aligned <- cor(
  z$mouse_direction_consistency,
  z$human_aligned_effect,
  method="spearman"
)
sp_i2_aligned <- cor(
  z$mouse_I2_percent,
  z$human_aligned_effect,
  method="spearman"
)
sp_mouse_human_effect <- cor(
  z$mouse_log2FC,
  z$human_RE_log2FC,
  method="spearman"
)
pe_mouse_human_effect <- cor(
  z$mouse_log2FC,
  z$human_RE_log2FC,
  method="pearson"
)

# PI descriptive comparison if available.
pi_available <- any(is.finite(z$mouse_PI_lower) & is.finite(z$mouse_PI_upper))
if (pi_available) {
  pi_tab <- z[, .(
    N=.N,
    Human_same_direction_fraction=mean(concordant_primary),
    Median_human_aligned_effect=median(human_aligned_effect)
  ), by=mouse_PI_excludes_zero]
} else {
  pi_tab <- data.table()
}

# Stable core descriptive only.
coretab <- z[mouse_stable_core==TRUE, .(
  mouse_ensembl, mouse_symbol,
  human_ENTREZID, human_symbol,
  mouse_log2FC, mouse_nGSE,
  mouse_direction_consistency, mouse_I2_percent,
  human_RE_log2FC, human_RE_FDR,
  concordant_primary,
  human_aligned_effect
)]
fwrite(coretab, COREOUT, sep="\t", quote=FALSE, na="NA")

# Sort final joined table without data.table inline expression bug.
z[, abs_mouse_effect_sort := abs(mouse_log2FC)]
setorder(z, -mouse_direction_consistency, abs_mouse_effect_sort)
z[, abs_mouse_effect_sort := NULL]
fwrite(z, JOINED, sep="\t", quote=FALSE, na="NA", compress="gzip")

primary_model <- models[Model=="PRIMARY_unadjusted_RE_same_direction"][1]
adjusted_model <- models[Model=="PRIMARY_adjusted_for_mouse_effect_and_coverage"][1]
cov20_model <- models[Model=="SENSITIVITY_mouse_coverage_ge20_GSE"][1]
fe_model <- models[Model=="SENSITIVITY_human_FE_direction"][1]
strict_model <- models[Model=="SENSITIVITY_GSE130970_strict4_human_RE_direction"][1]

summary <- c(
  "===== MOUSE HFD -> HUMAN MASLD REPRODUCIBILITY CONCORDANCE SUMMARY v1.6 =====",
  "",
  "STATUS=PASS",
  "SUPERSEDES=v1_to_v1_3_service_retrieval_failures_v1_4_summary_bug_v1_5_join_key_type_bug",
  "FROZEN_MOUSE_META_INPUT=gene_level_v1_2",
  "FROZEN_MOUSE_STABLE_CORE_INPUT=stable_core_v1",
  "FROZEN_HUMAN_META_INPUT=gene_meta_v1",
  paste0("R_version=",R.version.string),
    paste0("org.Hs.eg.db_version=",as.character(packageVersion("org.Hs.eg.db"))),
  "ENSEMBL_RELEASE=116",
  "ENSEMBL_RELEASE_DATE=June_2026",
  "ENSEMBL_ORTHOLOG_SOURCE=v1_4_frozen_map_from_official_release_116_Compara_FTP",
  "ORTHOLOG_NETWORK_ACCESS_IN_v1_6=NO",
  "JOIN_KEY_TYPES=character_for_mouse_ensembl_human_ensembl_human_ENTREZID",
  paste0("ORTHOLOG_SOURCE_FILE=",sub(paste0("^",ROOT,"/"),"",ORTHO_V1_4)),
  paste0("ORTHOLOG_SOURCE_SHA256=",orth_source_sha256),
  paste0("ORTHOLOG_v1_5_SHA256=",orth_v1_5_sha256),
  "ORTHOLOG_RULE=high_confidence_ortholog_one2one",
  "TRANSLATIONAL_UNIVERSE=mouse_protein_coding_meta_eligible_x_human_3cohort_meta_eligible_x_one2one_ortholog",
  "",
  paste0("Mouse_protein_coding_meta_genes=",nrow(m)),
  paste0("Frozen_highconfidence_one2one_orthologs=",nrow(orth)),
  paste0("Final_joined_translational_genes=",nrow(z)),
  paste0("Final_genes_mouse_nGSE_ge20=",sum(z$mouse_nGSE>=20)),
  "",
  "PRIMARY_ENDPOINT=human_RE_effect_same_direction_as_mouse_pooled_effect",
  "PRIMARY_PREDICTOR=mouse_GSE_balanced_direction_consistency_per_10_percentage_points",
  "PRIMARY_MODEL=binomial_logistic_regression",
  "PRIMARY_SIGNIFICANCE_TEST=unadjusted_primary_model",
  "ADJUSTED_SUPPORTIVE_MODEL=plus_abs_mouse_log2FC_plus_log10_mouse_nGSE",
  "SECONDARY_PREDICTORS_AND_ENDPOINTS=exploratory_or_supportive",
  "",
  paste0("Overall_human_same_direction_fraction=",format(mean(z$concordant_primary),digits=6)),
  paste0("Overall_human_FDR05_same_direction_fraction=",format(mean(z$human_significant_concordant),digits=6)),
  paste0("Overall_all3_human_concordant_fraction=",format(mean(z$human_all3_concordant),digits=6)),
  "",
  paste0("Primary_OR_per_10pp=",format(primary_model$OR,digits=6)),
  paste0("Primary_95CI=",format(primary_model$OR_95CI_low,digits=6),"_to_",format(primary_model$OR_95CI_high,digits=6)),
  paste0("Primary_P=",format(primary_model$P,digits=6)),
  paste0("Adjusted_OR_per_10pp=",format(adjusted_model$OR,digits=6)),
  paste0("Adjusted_P=",format(adjusted_model$P,digits=6)),
  paste0("Coverage_ge20_OR_per_10pp=",format(cov20_model$OR,digits=6)),
  paste0("Coverage_ge20_P=",format(cov20_model$P,digits=6)),
  paste0("Human_FE_sensitivity_OR_per_10pp=",format(fe_model$OR,digits=6)),
  paste0("Human_FE_sensitivity_P=",format(fe_model$P,digits=6)),
  paste0("Strict4_sensitivity_OR_per_10pp=",format(strict_model$OR,digits=6)),
  paste0("Strict4_sensitivity_P=",format(strict_model$P,digits=6)),
  "",
  paste0("Spearman_mouse_direction_consistency_vs_human_aligned_effect=",format(sp_dir_aligned,digits=6)),
  paste0("Spearman_mouse_I2_vs_human_aligned_effect=",format(sp_i2_aligned,digits=6)),
  paste0("Mouse_vs_human_effect_Pearson_r=",format(pe_mouse_human_effect,digits=6)),
  paste0("Mouse_vs_human_effect_Spearman_rho=",format(sp_mouse_human_effect,digits=6)),
  "",
  "REPRODUCIBILITY_QUARTILES:"
)
for (i in seq_len(nrow(quart))) {
  rr <- quart[i]
  summary <- c(summary,paste0(
    rr$reproducibility_quartile,
    "|N=",rr$N,
    "|median_mouse_consistency=",format(rr$Median_mouse_direction_consistency,digits=4),
    "|human_same_direction=",format(rr$Human_same_direction_fraction,digits=5),
    "|human_FDR05_same_direction=",format(rr$Human_FDR05_same_direction_fraction,digits=5),
    "|median_aligned_human_effect=",format(rr$Median_human_aligned_effect,digits=5)
  ))
}

summary <- c(
  summary,"",
  paste0("Mouse_PI_available=",ifelse(pi_available,"YES","NO")),
  paste0("Stable_core_orthologs_in_translational_universe=",nrow(coretab)),
  "STABLE_CORE_ANALYSIS=DESCRIPTIVE_ONLY_due_small_protein_coding_subset",
  "",
  paste0("Ortholog_map=",sub(paste0("^",ROOT,"/"),"",ORTHO)),
  paste0("Joined_gene_table=",sub(paste0("^",ROOT,"/"),"",JOINED)),
  paste0("Quartile_table=",sub(paste0("^",ROOT,"/"),"",QUART)),
  paste0("Model_table=",sub(paste0("^",ROOT,"/"),"",MODELS)),
  paste0("Stable_core_table=",sub(paste0("^",ROOT,"/"),"",COREOUT)),
  "",
  "CAUSAL_INTERPRETATION=NO",
  "HUMAN_SIGNIFICANCE_REQUIRED_FOR_PRIMARY_CONCORDANCE=NO",
  "RESULT_DEPENDENT_ORTHOLOG_SELECTION=NO"
)

writeLines(summary,SUMMARY)
cat("\n",paste(summary,collapse="\n"),"\n",sep="")

