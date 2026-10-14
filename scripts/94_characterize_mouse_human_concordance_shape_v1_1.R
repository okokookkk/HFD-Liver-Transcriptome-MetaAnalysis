#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required packages: ", paste(miss, collapse=","))
})
library(data.table)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

IN <- file.path(
  ROOT, "07_human_validation/mouse_human_concordance_v1_6",
  "mouse_HFD_human_MASLD_one2one_ortholog_concordance_v1_6.tsv.gz"
)

OUTDIR <- file.path(
  ROOT, "07_human_validation/mouse_human_concordance_followup_v1_1"
)
METADIR <- file.path(
  ROOT,
  "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_followup_v1_1"
)
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(METADIR, recursive=TRUE, showWarnings=FALSE)

OUT_Q <- file.path(
  OUTDIR, "mouse_reproducibility_quartile_contrast_models_v1_1.tsv"
)
OUT_TOP <- file.path(
  OUTDIR, "mouse_top_reproducibility_enrichment_v1_1.tsv"
)
OUT_PI <- file.path(
  OUTDIR, "mouse_prediction_interval_human_concordance_v1_1.tsv"
)
OUT_SUM <- file.path(
  METADIR, "mouse_human_concordance_shape_summary_v1_1.txt"
)

if (!file.exists(IN)) stop("Missing frozen input: ", IN)

z <- fread(IN)

need <- c(
  "mouse_direction_consistency",
  "mouse_abs_log2FC",
  "mouse_nGSE",
  "mouse_PI_excludes_zero",
  "concordant_primary",
  "human_significant_concordant",
  "human_all3_concordant",
  "human_aligned_effect"
)
miss <- setdiff(need, names(z))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse=","))

z[, concordant_primary := as.logical(concordant_primary)]
z[, human_significant_concordant := as.logical(human_significant_concordant)]
z[, human_all3_concordant := as.logical(human_all3_concordant)]
z[, mouse_PI_excludes_zero := as.logical(mouse_PI_excludes_zero)]

if (max(z$mouse_direction_consistency, na.rm=TRUE) > 1.5)
  z[, mouse_direction_consistency := mouse_direction_consistency/100]

z[, log10_mouse_nGSE := log10(mouse_nGSE)]

# -------------------------------------------------------------
# Deterministic rank-balanced quartiles, same construction as 93.
# -------------------------------------------------------------
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

glm_extract <- function(fit, term, label) {
  sm <- summary(fit)$coefficients
  if (!(term %in% rownames(sm))) stop(label, ": term not found: ", term)
  b <- sm[term,"Estimate"]
  se <- sm[term,"Std. Error"]
  data.table(
    Model=label,
    Term=term,
    N=nobs(fit),
    Beta_logOR=b,
    SE=se,
    Z=sm[term,"z value"],
    P=sm[term,"Pr(>|z|)"],
    OR=exp(b),
    OR_95CI_low=exp(b-1.96*se),
    OR_95CI_high=exp(b+1.96*se)
  )
}

# -------------------------------------------------------------
# 1) Quartile contrasts, Q1 as reference.
# These are characterization only; primary continuous model remains 93.
# -------------------------------------------------------------
fit_q <- glm(
  concordant_primary ~ reproducibility_quartile,
  family=binomial(),
  data=z
)
fit_q_adj <- glm(
  concordant_primary ~ reproducibility_quartile +
    mouse_abs_log2FC + log10_mouse_nGSE,
  family=binomial(),
  data=z
)

qmods <- rbindlist(list(
  glm_extract(
    fit_q,
    "reproducibility_quartileQ2",
    "UNADJUSTED_Q2_vs_Q1"
  ),
  glm_extract(
    fit_q,
    "reproducibility_quartileQ3",
    "UNADJUSTED_Q3_vs_Q1"
  ),
  glm_extract(
    fit_q,
    "reproducibility_quartileQ4_high",
    "UNADJUSTED_Q4_vs_Q1"
  ),
  glm_extract(
    fit_q_adj,
    "reproducibility_quartileQ2",
    "ADJUSTED_Q2_vs_Q1"
  ),
  glm_extract(
    fit_q_adj,
    "reproducibility_quartileQ3",
    "ADJUSTED_Q3_vs_Q1"
  ),
  glm_extract(
    fit_q_adj,
    "reproducibility_quartileQ4_high",
    "ADJUSTED_Q4_vs_Q1"
  )
))
qmods[, FDR_BH := p.adjust(P, method="BH")]
fwrite(qmods, OUT_Q, sep="\t", quote=FALSE, na="NA")

# -------------------------------------------------------------
# 2) Highest reproducibility quartile vs lower 75%.
# Characterizes whether association is concentrated at the top.
# -------------------------------------------------------------
z[, top_Q4 := reproducibility_quartile == "Q4_high"]

fit_top <- glm(
  concordant_primary ~ top_Q4,
  family=binomial(),
  data=z
)
fit_top_adj <- glm(
  concordant_primary ~ top_Q4 +
    mouse_abs_log2FC + log10_mouse_nGSE,
  family=binomial(),
  data=z
)

topmods <- rbindlist(list(
  glm_extract(fit_top, "top_Q4TRUE", "Q4_vs_lower75_unadjusted"),
  glm_extract(fit_top_adj, "top_Q4TRUE", "Q4_vs_lower75_adjusted")
))
fwrite(topmods, OUT_TOP, sep="\t", quote=FALSE, na="NA")

# -------------------------------------------------------------
# 3) Prediction-interval exclusion as prespecified supportive feature.
# This is secondary/descriptive because PI-excluding genes are sparse.
# -------------------------------------------------------------
pi <- z[!is.na(mouse_PI_excludes_zero)]

if (nrow(pi) == 0L)
  stop("No genes have non-missing mouse PI exclusion status")
if (uniqueN(pi$mouse_PI_excludes_zero) != 2L)
  stop("Mouse PI exclusion status does not contain both TRUE and FALSE classes")

pi_summary <- pi[, .(
  N=.N,
  Human_same_direction_fraction=mean(concordant_primary),
  Human_FDR05_same_direction_fraction=mean(human_significant_concordant),
  Human_all3_concordant_fraction=mean(human_all3_concordant),
  Median_human_aligned_effect=median(human_aligned_effect)
), by=mouse_PI_excludes_zero]

fit_pi <- glm(
  concordant_primary ~ mouse_PI_excludes_zero,
  family=binomial(),
  data=pi
)
fit_pi_adj <- glm(
  concordant_primary ~ mouse_PI_excludes_zero +
    mouse_abs_log2FC + log10_mouse_nGSE,
  family=binomial(),
  data=pi
)

pimods <- rbindlist(list(
  glm_extract(
    fit_pi,
    "mouse_PI_excludes_zeroTRUE",
    "PI_excludes_zero_unadjusted"
  ),
  glm_extract(
    fit_pi_adj,
    "mouse_PI_excludes_zeroTRUE",
    "PI_excludes_zero_adjusted"
  )
))

fwrite(pi_summary, OUT_PI, sep="\t", quote=FALSE, na="NA")
fwrite(
  pimods,
  sub("\\.tsv$","_models.tsv",OUT_PI),
  sep="\t", quote=FALSE, na="NA"
)

# -------------------------------------------------------------
# Summary.
# -------------------------------------------------------------
qsum <- z[, .(
  N=.N,
  Median_mouse_consistency=median(mouse_direction_consistency),
  Human_same_direction_fraction=mean(concordant_primary),
  Human_FDR05_same_direction_fraction=mean(human_significant_concordant),
  Median_human_aligned_effect=median(human_aligned_effect)
), by=reproducibility_quartile]
setorder(qsum, reproducibility_quartile)

q4u <- qmods[Model=="UNADJUSTED_Q4_vs_Q1"][1]
q4a <- qmods[Model=="ADJUSTED_Q4_vs_Q1"][1]
topu <- topmods[Model=="Q4_vs_lower75_unadjusted"][1]
topa <- topmods[Model=="Q4_vs_lower75_adjusted"][1]
piu <- pimods[Model=="PI_excludes_zero_unadjusted"][1]
pia <- pimods[Model=="PI_excludes_zero_adjusted"][1]

summary <- c(
  "===== MOUSE-HUMAN CONCORDANCE SHAPE SUMMARY v1.1 =====",
  "",
  "STATUS=PASS",
  "SUPERSEDES=v1_PI_output_merge_bug_only",
  "UPSTREAM_PRIMARY_ANALYSIS=93_v1_6",
  "PRIMARY_INFERENCE_CHANGED=NO",
  "FOLLOWUP_ROLE=secondary_characterization_of_effect_shape",
  "",
  paste0("Joined_genes=",nrow(z)),
  "",
  "QUARTILE_DESCRIPTIVES:"
)

for (i in seq_len(nrow(qsum))) {
  rr <- qsum[i]
  summary <- c(summary, paste0(
    rr$reproducibility_quartile,
    "|N=",rr$N,
    "|median_mouse_consistency=",format(rr$Median_mouse_consistency,digits=5),
    "|human_same_direction=",format(rr$Human_same_direction_fraction,digits=6),
    "|human_FDR05_same_direction=",format(rr$Human_FDR05_same_direction_fraction,digits=6),
    "|median_human_aligned_effect=",format(rr$Median_human_aligned_effect,digits=6)
  ))
}

summary <- c(
  summary,
  "",
  paste0("Q4_vs_Q1_unadjusted_OR=",format(q4u$OR,digits=6)),
  paste0(
    "Q4_vs_Q1_unadjusted_95CI=",
    format(q4u$OR_95CI_low,digits=6),"_to_",
    format(q4u$OR_95CI_high,digits=6)
  ),
  paste0("Q4_vs_Q1_unadjusted_P=",format(q4u$P,digits=6)),
  paste0("Q4_vs_Q1_adjusted_OR=",format(q4a$OR,digits=6)),
  paste0("Q4_vs_Q1_adjusted_P=",format(q4a$P,digits=6)),
  "",
  paste0("Q4_vs_lower75_unadjusted_OR=",format(topu$OR,digits=6)),
  paste0("Q4_vs_lower75_unadjusted_P=",format(topu$P,digits=6)),
  paste0("Q4_vs_lower75_adjusted_OR=",format(topa$OR,digits=6)),
  paste0("Q4_vs_lower75_adjusted_P=",format(topa$P,digits=6)),
  "",
  paste0("PI_excludes_zero_N=",sum(pi$mouse_PI_excludes_zero)),
  paste0("PI_crosses_zero_N=",sum(!pi$mouse_PI_excludes_zero)),
  paste0("PI_excludes_zero_unadjusted_OR=",format(piu$OR,digits=6)),
  paste0("PI_excludes_zero_unadjusted_P=",format(piu$P,digits=6)),
  paste0("PI_excludes_zero_adjusted_OR=",format(pia$OR,digits=6)),
  paste0("PI_excludes_zero_adjusted_P=",format(pia$P,digits=6)),
  "",
  "INTERPRETATION_RULE=do_not_replace_prespecified_continuous_primary_model_with_posthoc_threshold",
  "CAUSAL_INTERPRETATION=NO"
)

writeLines(summary, OUT_SUM)
cat("\n",paste(summary,collapse="\n"),"\n",sep="")

