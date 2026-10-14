#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Required R package not installed: DESeq2")
  }
})

library(DESeq2)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

AUDIT <- file.path(
  ROOT,
  "00_metadata/analysis_inputs/frozen_v1",
  "HFD_primary105_contrast_analysis_input_audit_v1_1.tsv"
)

FINAL_AUDIT <- file.path(
  ROOT,
  "00_metadata/processing_production/final_audit_v1",
  "HFD_primary64_final_production_audit_summary_v1.txt"
)

OUTROOT <- file.path(ROOT, "04_de", "contrast_level_v1")
METADIR <- file.path(
  ROOT,
  "00_metadata/analysis_outputs/contrast_de_v1"
)

dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

STATUS_OUT <- file.path(METADIR, "HFD_primary105_DESeq2_status_v1.tsv")
INDEX_OUT <- file.path(METADIR, "HFD_primary105_DESeq2_effect_file_index_v1.tsv")
SUMMARY_OUT <- file.path(METADIR, "HFD_primary105_DESeq2_summary_v1.txt")

if (!file.exists(AUDIT)) stop("Missing analysis-input audit: ", AUDIT)
if (!file.exists(FINAL_AUDIT)) stop("Missing final production audit: ", FINAL_AUDIT)

audit <- read.delim(
  AUDIT,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_cols <- c(
  "Contrast_id", "GSE", "N_HFD_GSM", "N_CONTROL_GSM",
  "HFD_GSMs", "CONTROL_GSMs", "Final_matrix", "Constructable"
)

missing_cols <- setdiff(required_cols, colnames(audit))
if (length(missing_cols) > 0) {
  stop("Analysis-input audit missing columns: ", paste(missing_cols, collapse = ", "))
}

if (nrow(audit) != 105L || length(unique(audit$Contrast_id)) != 105L) {
  stop("Expected exactly 105 unique contrasts")
}

if (!all(audit$Constructable == "PASS")) {
  stop("At least one frozen contrast is not constructable")
}

split_gsms <- function(s) {
  z <- trimws(unlist(strsplit(s, ",", fixed = TRUE)))
  z[nzchar(z)]
}

atomic_write_tsv_gz <- function(df, outfile) {
  tmp <- paste0(outfile, ".tmp")
  con <- gzfile(tmp, open = "wt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  close(con)
  on.exit(NULL, add = FALSE)

  if (file.exists(outfile)) file.remove(outfile)
  if (!file.rename(tmp, outfile)) {
    stop("Failed atomic rename: ", tmp, " -> ", outfile)
  }
}

status_rows <- list()
index_rows <- list()
n_fail <- 0L

cat("===== HFD PRIMARY105 CONTRAST-LEVEL DESeq2 v1 =====\n")
cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n\n")

for (ii in seq_len(nrow(audit))) {

  r <- audit[ii, , drop = FALSE]
  cid <- r$Contrast_id
  gse <- r$GSE

  cat(sprintf("[%03d/105] %s (%s)\n", ii, cid, gse))

  hfd <- split_gsms(r$HFD_GSMs)
  ctl <- split_gsms(r$CONTROL_GSMs)

  contrast_outdir <- file.path(OUTROOT, gse)
  dir.create(contrast_outdir, recursive = TRUE, showWarnings = FALSE)

  effect_file <- file.path(
    contrast_outdir,
    paste0(cid, "_DESeq2_effects_v1.tsv.gz")
  )
  done_file <- file.path(
    contrast_outdir,
    paste0(cid, "_DESeq2_DONE_v1.txt")
  )

  # Restart-safe behavior: completed versioned outputs are never overwritten.
  if (file.exists(effect_file) && file.exists(done_file)) {
    dl <- readLines(done_file, warn = FALSE)
    kv <- strsplit(dl[grepl("=", dl, fixed = TRUE)], "=", fixed = TRUE)
    keys <- vapply(kv, `[`, character(1), 1L)
    vals <- vapply(kv, function(z) paste(z[-1L], collapse = "="), character(1))
    done_map <- setNames(vals, keys)

    if (!identical(unname(done_map["STATUS"]), "PASS")) {
      stop(cid, ": existing DONE marker is not PASS")
    }

    get_int <- function(k) {
      v <- unname(done_map[k])
      if (length(v) == 0L || is.na(v) || !nzchar(v)) return(NA_integer_)
      as.integer(v)
    }

    cat("  SKIP: existing PASS effect file + DONE marker\n")

    status_rows[[length(status_rows) + 1L]] <- data.frame(
      Contrast_id = cid,
      GSE = gse,
      Status = "PASS",
      N_CONTROL_GSM = length(ctl),
      N_HFD_GSM = length(hfd),
      Input_gene_rows = get_int("Input_gene_rows"),
      Modeled_genes = get_int("Modeled_genes"),
      Finite_effect_SE = get_int("Finite_effect_SE"),
      Cooks_flagged_effects = get_int("Cooks_flagged_effects"),
      Primary_effect_QC_PASS = get_int("Primary_effect_QC_PASS"),
      Error = "",
      stringsAsFactors = FALSE
    )

    index_rows[[length(index_rows) + 1L]] <- data.frame(
      Contrast_id = cid,
      GSE = gse,
      Effect_file = sub(paste0("^", ROOT, "/"), "", effect_file),
      DONE_file = sub(paste0("^", ROOT, "/"), "", done_file),
      stringsAsFactors = FALSE
    )

    next
  }

  one <- tryCatch({

    if (length(hfd) != as.integer(r$N_HFD_GSM)) {
      stop("HFD GSM count mismatch")
    }
    if (length(ctl) != as.integer(r$N_CONTROL_GSM)) {
      stop("CONTROL GSM count mismatch")
    }
    if (length(intersect(hfd, ctl)) > 0) {
      stop("HFD and CONTROL GSM sets overlap")
    }
    if (length(hfd) < 3L || length(ctl) < 3L) {
      stop("Frozen contrast unexpectedly has fewer than 3 GSM in an arm")
    }

    matrix_file <- file.path(ROOT, r$Final_matrix)
    if (!file.exists(matrix_file)) {
      stop("Final GSM matrix missing: ", matrix_file)
    }

    x <- read.delim(
      matrix_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = ""
    )

    if (!"Geneid" %in% colnames(x)) stop("Final matrix lacks Geneid")
    if (anyDuplicated(x$Geneid)) stop("Duplicate Geneid in final matrix")

    samples <- c(ctl, hfd)
    missing_samples <- setdiff(samples, colnames(x))
    if (length(missing_samples) > 0) {
      stop(
        "Frozen GSM missing from final matrix: ",
        paste(missing_samples, collapse = ",")
      )
    }

    count_df <- x[, samples, drop = FALSE]

    # Validate nonnegative integer-like raw counts.
    count_mat_num <- as.matrix(count_df)
    storage.mode(count_mat_num) <- "numeric"

    if (any(!is.finite(count_mat_num))) stop("Non-finite raw counts")
    if (any(count_mat_num < 0)) stop("Negative raw counts")
    if (any(abs(count_mat_num - round(count_mat_num)) > 1e-8)) {
      stop("Non-integer raw counts")
    }
    if (any(count_mat_num > .Machine$integer.max)) {
      stop("Raw count exceeds R integer range")
    }

    count_mat <- matrix(
      as.integer(round(count_mat_num)),
      nrow = nrow(count_mat_num),
      ncol = ncol(count_mat_num),
      dimnames = list(x$Geneid, samples)
    )

    # Frozen SAP v1 expression filter:
    # >=10 raw counts in >=3 biological samples within this contrast.
    n_ge10 <- rowSums(count_mat >= 10L)
    keep <- n_ge10 >= 3L

    n_input_genes <- nrow(count_mat)
    n_modeled_genes <- sum(keep)

    if (n_modeled_genes < 100L) {
      stop("Unexpectedly few genes after expression filtering: ", n_modeled_genes)
    }

    count_keep <- count_mat[keep, , drop = FALSE]
    n_ge10_keep <- n_ge10[keep]

    group <- factor(
      c(rep("CONTROL", length(ctl)), rep("HFD", length(hfd))),
      levels = c("CONTROL", "HFD")
    )

    coldata <- data.frame(
      row.names = samples,
      group = group
    )

    dds <- DESeqDataSetFromMatrix(
      countData = count_keep,
      colData = coldata,
      design = ~ group
    )

    # Disable automatic count replacement. We retain Cook's-distance evidence
    # and perform explicit effect-level QC instead of silently replacing counts.
    dds <- DESeq(
      dds,
      test = "Wald",
      betaPrior = FALSE,
      minReplicatesForReplace = Inf,
      quiet = TRUE
    )

    # Primary meta-analysis uses unshrunken MLE log2FC + Wald SE.
    # Independent filtering is disabled because within-contrast significance
    # does not determine meta-analysis eligibility.
    res <- results(
      dds,
      contrast = c("group", "HFD", "CONTROL"),
      independentFiltering = FALSE,
      cooksCutoff = FALSE,
      alpha = 0.05
    )

    res_df <- as.data.frame(res)
    res_df$Geneid <- rownames(res_df)

    # Group-specific normalized means are descriptive QC variables only.
    norm <- counts(dds, normalized = TRUE)
    mean_ctl <- rowMeans(norm[, ctl, drop = FALSE])
    mean_hfd <- rowMeans(norm[, hfd, drop = FALSE])

    # Explicitly reproduce the standard DESeq2 Cook's-distance cutoff so that
    # contrast-gene effects can be flagged without using within-study p-values.
    cooks <- assays(dds)[["cooks"]]
    if (is.null(cooks)) {
      max_cooks <- rep(NA_real_, nrow(dds))
      names(max_cooks) <- rownames(dds)
      cooks_cutoff <- NA_real_
    } else {
      max_cooks <- apply(cooks, 1, function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0L) return(NA_real_)
        max(v)
      })
      mm <- model.matrix(design(dds), colData(dds))
      p_model <- ncol(mm)
      m_samples <- nrow(colData(dds))
      cooks_cutoff <- qf(0.99, p_model, m_samples - p_model)
    }

    out <- data.frame(
      Geneid = res_df$Geneid,
      baseMean = res_df$baseMean,
      mean_normalized_CONTROL = mean_ctl[res_df$Geneid],
      mean_normalized_HFD = mean_hfd[res_df$Geneid],
      log2FoldChange = res_df$log2FoldChange,
      lfcSE = res_df$lfcSE,
      sampling_variance = res_df$lfcSE^2,
      Wald_stat = res_df$stat,
      Wald_pvalue = res_df$pvalue,
      N_samples_count_ge_10 = n_ge10_keep[res_df$Geneid],
      N_CONTROL_GSM = length(ctl),
      N_HFD_GSM = length(hfd),
      max_Cooks_distance = max_cooks[res_df$Geneid],
      DESeq2_Cooks_cutoff = cooks_cutoff,
      DESeq2_Cooks_flag = ifelse(
        is.finite(max_cooks[res_df$Geneid]) &
          is.finite(cooks_cutoff) &
          max_cooks[res_df$Geneid] > cooks_cutoff,
        "YES", "NO"
      ),
      stringsAsFactors = FALSE
    )

    out$Finite_effect_SE <- ifelse(
      is.finite(out$log2FoldChange) &
        is.finite(out$lfcSE) &
        out$lfcSE > 0,
      "YES", "NO"
    )

    # Primary effect-level QC for downstream meta-analysis:
    # finite effect/SE and no standard DESeq2 Cook's-distance outlier flag.
    out$Primary_effect_QC <- ifelse(
      out$Finite_effect_SE == "YES" &
        out$DESeq2_Cooks_flag == "NO",
      "PASS", "FAIL"
    )

    # Stable sort by the common gene order in the final GSM matrix.
    ord <- match(out$Geneid, x$Geneid)
    out <- out[order(ord), , drop = FALSE]

    atomic_write_tsv_gz(out, effect_file)

    n_finite <- sum(out$Finite_effect_SE == "YES")
    n_cooks <- sum(out$DESeq2_Cooks_flag == "YES", na.rm = TRUE)
    n_primary_pass <- sum(out$Primary_effect_QC == "PASS")

    done_lines <- c(
      paste0("Contrast_id=", cid),
      paste0("GSE=", gse),
      "STATUS=PASS",
      paste0("N_CONTROL_GSM=", length(ctl)),
      paste0("N_HFD_GSM=", length(hfd)),
      paste0("Input_gene_rows=", n_input_genes),
      paste0("Modeled_genes=", n_modeled_genes),
      paste0("Finite_effect_SE=", n_finite),
      paste0("Cooks_flagged_effects=", n_cooks),
      paste0("Primary_effect_QC_PASS=", n_primary_pass),
      "Expression_filter=raw_count_ge_10_in_ge_3_biological_samples",
      "Effect_size=unshrunken_DESeq2_log2FoldChange",
      "Sampling_variance=lfcSE_squared",
      "Independent_filtering=FALSE",
      "Automatic_count_replacement=DISABLED",
      "Meta_analysis_performed=NO",
      paste0("Effect_file=", sub(paste0("^", ROOT, "/"), "", effect_file))
    )
    writeLines(done_lines, con = done_file)

    list(
      ok = TRUE,
      n_input = n_input_genes,
      n_modeled = n_modeled_genes,
      n_finite = n_finite,
      n_cooks = n_cooks,
      n_pass = n_primary_pass,
      error = ""
    )

  }, error = function(e) {
    list(
      ok = FALSE,
      n_input = NA_integer_,
      n_modeled = NA_integer_,
      n_finite = NA_integer_,
      n_cooks = NA_integer_,
      n_pass = NA_integer_,
      error = conditionMessage(e)
    )
  })

  if (!one$ok) {
    n_fail <- n_fail + 1L
    cat("  FAIL:", one$error, "\n")
  } else {
    cat(
      sprintf(
        "  PASS: modeled=%d finite=%d cooks_flagged=%d primary_QC=%d\n",
        one$n_modeled, one$n_finite, one$n_cooks, one$n_pass
      )
    )
  }

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    Contrast_id = cid,
    GSE = gse,
    Status = ifelse(one$ok, "PASS", "FAIL"),
    N_CONTROL_GSM = length(ctl),
    N_HFD_GSM = length(hfd),
    Input_gene_rows = one$n_input,
    Modeled_genes = one$n_modeled,
    Finite_effect_SE = one$n_finite,
    Cooks_flagged_effects = one$n_cooks,
    Primary_effect_QC_PASS = one$n_pass,
    Error = one$error,
    stringsAsFactors = FALSE
  )

  index_rows[[length(index_rows) + 1L]] <- data.frame(
    Contrast_id = cid,
    GSE = gse,
    Effect_file = sub(paste0("^", ROOT, "/"), "", effect_file),
    DONE_file = sub(paste0("^", ROOT, "/"), "", done_file),
    stringsAsFactors = FALSE
  )
}

status <- do.call(rbind, status_rows)
index <- do.call(rbind, index_rows)

write.table(
  status,
  STATUS_OUT,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  index,
  INDEX_OUT,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

n_pass_contrasts <- sum(status$Status == "PASS")
n_fail_contrasts <- sum(status$Status == "FAIL")

summary <- c(
  "===== HFD PRIMARY105 CONTRAST-LEVEL DESeq2 SUMMARY v1 =====",
  "",
  paste0("FINAL_STATUS=", ifelse(n_fail_contrasts == 0L, "PASS", "FAIL")),
  paste0("Contrasts_total=", nrow(status)),
  paste0("Contrasts_PASS=", n_pass_contrasts),
  paste0("Contrasts_FAIL=", n_fail_contrasts),
  paste0("Independent_GSE=", length(unique(status$GSE))),
  "",
  paste0("R_version=", R.version.string),
  paste0("DESeq2_version=", as.character(packageVersion("DESeq2"))),
  "",
  paste0("Modeled_genes_min=", min(status$Modeled_genes, na.rm = TRUE)),
  paste0("Modeled_genes_median=", median(status$Modeled_genes, na.rm = TRUE)),
  paste0("Modeled_genes_max=", max(status$Modeled_genes, na.rm = TRUE)),
  paste0("Finite_effect_SE_min=", min(status$Finite_effect_SE, na.rm = TRUE)),
  paste0("Finite_effect_SE_median=", median(status$Finite_effect_SE, na.rm = TRUE)),
  paste0("Finite_effect_SE_max=", max(status$Finite_effect_SE, na.rm = TRUE)),
  paste0("Cooks_flagged_effects_total=", sum(status$Cooks_flagged_effects, na.rm = TRUE)),
  paste0("Primary_effect_QC_PASS_total=", sum(status$Primary_effect_QC_PASS, na.rm = TRUE)),
  "",
  "Expression_filter=raw_count_ge_10_in_ge_3_biological_samples_per_contrast",
  "Primary_effect=unshrunken_log2FC",
  "Sampling_variance=lfcSE_squared",
  "Within_contrast_DE_significance_used_for_meta_eligibility=NO",
  "Independent_filtering=FALSE",
  "Automatic_DESeq2_count_replacement=DISABLED",
  "Cook_distance_flag_recorded=YES",
  "Meta_analysis_performed=NO",
  "",
  paste0("Status_table=", sub(paste0("^", ROOT, "/"), "", STATUS_OUT)),
  paste0("Effect_file_index=", sub(paste0("^", ROOT, "/"), "", INDEX_OUT))
)

writeLines(summary, con = SUMMARY_OUT)

cat("\n", paste(summary, collapse = "\n"), "\n", sep = "")

if (n_fail > 0L) {
  quit(status = 1L)
}

