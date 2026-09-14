# R/run_pca.R
# Sample-level PCA per subcohort (exploratory QC + visualization source).
#
# Design notes:
#   - Identical matrix construction + missingness rule as run_de.R
#     (>20% drop, median-impute): QC and DE share preprocessing.
#   - Proteins are z-scored ACROSS SAMPLES before PCA, so PCs reflect
#     correlation structure, not raw NPX abundance (high-abundance
#     assays would otherwise dominate PC1). Direction matters:
#     scale() scales COLUMNS, so we feed it t(mat) whose columns
#     are proteins. t(scale(t(mat))) would scale per SAMPLE -- wrong.
#   - Samples are the observations: prcomp() on the standardized
#     samples x proteins matrix.
#   - Plate is retained in the artifact so the app can color by plate
#     (batch-effect QC: do samples cluster by Plate_ID?).
#
# Reads:  data/harmonized_npx.rds (via caller)
# Writes: outputs/<subcohort>/pca.rds
#           list(result = scores data.frame, params = list)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(here)
})

run_pca <- function(npx, subcohort_name, out_dir,
                    max_missing = 0.20) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df <- npx |> filter(subcohort == .env$subcohort_name)

  # --- 1. Protein x sample matrix (same recipe as run_de.R) ----------------
  mat <- df |>
    select(Assay, SampleID, NPX) |>
    pivot_wider(names_from = SampleID, values_from = NPX) |>
    column_to_rownames("Assay") |>
    as.matrix()

  # --- 2. Missingness: >20% drop, median-impute the rest --------------------
  frac_missing <- rowMeans(is.na(mat))
  keep      <- frac_missing <= max_missing
  n_dropped <- sum(!keep)
  mat       <- mat[keep, , drop = FALSE]

  n_imputed <- sum(is.na(mat))
  if (n_imputed > 0) {
    med_by_row <- apply(mat, 1, median, na.rm = TRUE)
    idx        <- which(is.na(mat), arr.ind = TRUE)
    mat[idx]   <- med_by_row[idx[, "row"]]
  }
  message(sprintf(
    "[%s] PCA matrix: %d proteins x %d samples (%d dropped >%d%% missing; %d cells median-imputed)",
    subcohort_name, nrow(mat), ncol(mat), n_dropped,
    round(max_missing * 100), n_imputed
  ))

  # --- 3. Standardize per protein, samples as observations -------------------
  # z: samples x proteins, each protein column z-scored across samples.
  z  <- scale(t(mat))
  pr <- prcomp(z)   # center = TRUE is a no-op post-scaling; scale. = FALSE

  var_explained <- summary(pr)$importance[2, ]   # proportion per PC
  message(sprintf(
    "[%s] PCA: PC1 %.1f%% / PC2 %.1f%% / PC3 %.1f%% of variance",
    subcohort_name,
    var_explained[1] * 100, var_explained[2] * 100, var_explained[3] * 100
  ))

  # --- 4. Scores + aligned sample metadata ------------------------------------
  # pr$x columns are already named PC1, PC2, ... (307 PCs for plasma;
  # the app/report only consume the first few).
  scores <- as.data.frame(pr$x) |>
    rownames_to_column("SampleID") |>
    left_join(
      df |> distinct(SampleID, Individual, Group, Sex, Age, Severity, Plate),
      by = "SampleID"
    )

  params <- list(
    subcohort           = subcohort_name,
    method              = "prcomp (SVD), samples as observations",
    scaling             = "per-protein z-score across samples",
    max_missing         = max_missing,
    n_samples           = nrow(z),
    n_proteins_retained = ncol(z),
    n_proteins_dropped  = n_dropped,
    n_cells_imputed     = n_imputed,
    var_explained       = var_explained,          # full vector, named PCs
    var_pc1_pc2_cum     = var_explained[1] + var_explained[2],
    finished_utc        = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  saveRDS(list(result = scores, params = params),
          file.path(out_dir, "pca.rds"))
  invisible(list(result = scores, params = params))
}
