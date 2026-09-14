# R/run_de.R
# Differential expression (COVID POSITIVE vs NEGATIVE) per subcohort with limma.
#
# Design decisions (from data diagnostics, 2026-09-14):
#   - Longitudinal structure: samples repeated within individuals
#     (plasma: 307 samples / 105 individuals). Handled with
#     duplicateCorrelation() when enough individuals are replicated;
#     otherwise plain lmFit with a logged note.
#   - Covariates: Sex + (binned) Age for subcohorts with enough n;
#     passed explicitly per subcohort, with NA/rank guards.
#   - Missingness is encoded as ABSENT ROWS in the source; after
#     pivot_wider these become NAs. Rule: drop proteins >20% missing,
#     median-impute the scattered remainder. Counts are logged.
#   - Contrast: POSITIVE - NEGATIVE (verified group levels), so
#     logFC > 0 = higher in COVID POSITIVE.
#
# Reads:  data/harmonized_npx.rds (via caller)
# Writes: outputs/<subcohort>/de_results.rds
#           list(result = data.frame, params = list)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(here)
  library(limma)
})

run_de <- function(npx, subcohort_name, out_dir,
                   covariates = character(0),
                   max_missing    = 0.20,
                   min_rep_indivs = 10) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df <- npx |> filter(subcohort == .env$subcohort_name)

  # --- 1. Protein x sample matrix ------------------------------------------
  mat <- df |>
    select(Assay, SampleID, NPX) |>
    pivot_wider(names_from = SampleID, values_from = NPX) |>
    column_to_rownames("Assay") |>
    as.matrix()

  # --- 2. Missingness: >20% drop, median-impute the rest --------------------
  frac_missing <- rowMeans(is.na(mat))
  keep        <- frac_missing <= max_missing
  n_dropped   <- sum(!keep)
  mat         <- mat[keep, , drop = FALSE]

  n_imputed <- sum(is.na(mat))
  if (n_imputed > 0) {
    med_by_row <- apply(mat, 1, median, na.rm = TRUE)
    idx        <- which(is.na(mat), arr.ind = TRUE)
    mat[idx]   <- med_by_row[idx[, "row"]]
  }
  message(sprintf(
    "[%s] DE matrix: %d proteins x %d samples (%d dropped >%d%% missing; %d cells median-imputed)",
    subcohort_name, nrow(mat), ncol(mat), n_dropped,
    round(max_missing * 100), n_imputed
  ))

  # --- 3. Sample metadata, aligned to matrix columns -------------------------
  # NOTE: select(any_of()) + distinct(), NOT distinct(all_of()) --
  # in distinct(), all_of() is evaluated as a data-mask expression
  # ("size N or 1" error), while in select() it is tidyselect.
  req_cov    <- covariates
  covariates <- intersect(req_cov, names(df))
  bad_cov    <- setdiff(req_cov, covariates)
  if (length(bad_cov) > 0) {
    message(sprintf("[%s] DE: covariate(s) not in data, ignored: %s",
                    subcohort_name, paste(bad_cov, collapse = ", ")))
  }

  meta <- df |>
    select(SampleID, Individual, Group, any_of(covariates)) |>
    distinct() |>
    filter(SampleID %in% colnames(mat)) |>
    arrange(match(SampleID, colnames(mat)))

  # Covariate hygiene: drop any covariate containing NA (would break
  # model.matrix); convert character covariates to factor for
  # deterministic level order.
  keep_cov <- covariates[vapply(covariates, function(v) {
    !any(is.na(meta[[v]]))
  }, logical(1))]
  dropped_cov <- setdiff(covariates, keep_cov)
  if (length(dropped_cov) > 0) {
    message(sprintf("[%s] DE: covariate(s) with NA dropped: %s",
                    subcohort_name, paste(dropped_cov, collapse = ", ")))
  }
  for (v in keep_cov) {
    if (is.character(meta[[v]])) meta[[v]] <- factor(meta[[v]])
  }

  # --- 4. Design matrix -------------------------------------------------------
  meta$Group <- factor(meta$Group, levels = c("POSITIVE", "NEGATIVE"))
  build_design <- function(covs) {
    rhs <- paste(c("0 + Group", covs), collapse = " + ")
    model.matrix(as.formula(paste("~", rhs)), data = meta)
  }
  design <- build_design(keep_cov)

  # Rank guard: a covariate level confounded with Group (e.g. an age bin
  # containing only POSITIVE samples) makes the design rank-deficient.
  # Fall back to Group-only rather than fit garbage.
  if (qr(design)$rank < ncol(design)) {
    message(sprintf(
      "[%s] DE: design rank-deficient with covariates; falling back to Group only",
      subcohort_name))
    keep_cov <- character(0)
    design   <- build_design(keep_cov)
  }

  # --- 5. Fit: blocked (duplicateCorrelation) when replication supports it ---
  n_rep     <- sum(table(meta$Individual) >= 2)
  use_block <- n_rep >= min_rep_indivs

  consensus <- NULL
  if (use_block) {
    corfit    <- duplicateCorrelation(mat, design, block = meta$Individual)
    consensus <- corfit$consensus
    fit <- lmFit(mat, design, block = meta$Individual, correlation = consensus)
    message(sprintf(
      "[%s] DE: blocked fit via duplicateCorrelation (%d/%d individuals with >=2 samples; consensus r = %.3f)",
      subcohort_name, n_rep, n_distinct(meta$Individual), consensus
    ))
  } else {
    fit <- lmFit(mat, design)
    message(sprintf(
      "[%s] DE: plain fit (only %d individuals with >=2 samples < threshold %d); within-person correlation not estimable",
      subcohort_name, n_rep, min_rep_indivs
    ))
  }

  # --- 6. Contrast + empirical Bayes -----------------------------------------
  # Manual contrast vector instead of makeContrasts(): covariate factor
  # levels like "Age(40,60]" produce design column names that are not
  # syntactically valid R names, which makeContrasts() rejects. A hand-
  # built vector only needs the two Group columns, so any covariate
  # naming is safe. Entries: +1 = POSITIVE, -1 = NEGATIVE.
  cm <- matrix(0, nrow = ncol(design),
               dimnames = list(colnames(design), "POSITIVE_vs_NEGATIVE"))
  cm["GroupPOSITIVE", 1] <-  1
  cm["GroupNEGATIVE", 1] <- -1
  fit2 <- eBayes(contrasts.fit(fit, cm))

  # --- 7. Results table --------------------------------------------------------
  res <- topTable(fit2, number = Inf, sort.by = "P") |>
    rownames_to_column("Assay") |>
    left_join(df |> distinct(Assay, UniProt, Panel), by = "Assay") |>
    select(Assay, UniProt, Panel, logFC, AveExpr, t, P.Value, adj.P.Val, B)

  params <- list(
    subcohort                = subcohort_name,
    contrast                 = "POSITIVE - NEGATIVE (logFC > 0 = higher in COVID POSITIVE)",
    method                   = "limma lmFit + eBayes, BH FDR",
    blocked                  = use_block,
    consensus_correlation    = if (is.null(consensus)) NA_real_ else consensus,
    covariates               = keep_cov,
    min_rep_indivs_threshold = min_rep_indivs,
    max_missing              = max_missing,
    n_samples                = ncol(mat),
    n_individuals            = n_distinct(meta$Individual),
    n_proteins_tested        = nrow(mat),
    n_proteins_dropped       = n_dropped,
    n_cells_imputed          = n_imputed,
    n_sig_fdr05              = sum(res$adj.P.Val < 0.05),
    finished_utc             = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  saveRDS(list(result = res, params = params),
          file.path(out_dir, "de_results.rds"))
  message(sprintf("[%s] DE: %d proteins tested, %d with adj.P < 0.05",
                  subcohort_name, nrow(res), params$n_sig_fdr05))
  invisible(list(result = res, params = params))
}
