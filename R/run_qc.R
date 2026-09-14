# R/run_qc.R
# QC metrics per subcohort -> self-describing artifact.
#
# result (named list of data frames; QC tab consumes by name):
#   $assay_missingness   -- per-assay observed/missing/frac, >20% flag
#   $sample_missingness -- per-sample observed/missing/frac, >20% flag
#   $outliers            -- PCA-distance outlier screen (see below)
#   $plate_summary       -- samples per Plate x Group
#
# Outlier screen: Euclidean distance to the sample's OWN group centroid in
# PC1-PC5 (projection from outputs/<subcohort>/pca.rds -- pipeline reuses
# the artifact, does not recompute PCA). Robust z = (dist - median)/MAD
# within group, one-sided (distance can only be too large), flag z > 3.
# Screening within group means the disease signal cannot masquerade as
# outlyingness.
#
# Missingness is encoded as ABSENT rows in the source; per-unit counts are
# expected - observed. A negative count would indicate duplicated rows --
# guarded with a hard stop (harmonization verified uniqueness, so silence
# of this guard is itself a check).
#
# QC flags are INFORMATIONAL: no samples are excluded anywhere in this
# project (no pre-specified exclusion rule; see README - Methods).
# LOD flags are not computable (LOD fields absent from the harmonized
# dataset) -- recorded in params rather than silently omitted.
#
# Reads:  data/harmonized_npx.rds (via caller),
#         outputs/<subcohort>/pca.rds (hard prerequisite)
# Writes: outputs/<subcohort>/qc.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})

run_qc <- function(npx, subcohort_name, out_dir,
                   outlier_z = 3, n_pcs = 5, max_missing = 0.20) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df <- npx |> filter(subcohort == .env$subcohort_name)

  n_samples <- n_distinct(df$SampleID)
  n_assays  <- n_distinct(df$Assay)

  # --- 1. Per-assay + per-sample missingness --------------------------------
  assay_miss <- df |>
    group_by(Assay, Panel) |>
    summarise(n_observed = n(), .groups = "drop") |>
    mutate(
      n_missing    = n_samples - n_observed,
      frac_missing = n_missing / n_samples,
      flag_missing = frac_missing > max_missing
    ) |>
    arrange(desc(frac_missing), Assay)

  sample_miss <- df |>
    group_by(SampleID, Individual, Group) |>
    summarise(n_observed = n(), .groups = "drop") |>
    mutate(
      n_missing    = n_assays - n_observed,
      frac_missing = n_missing / n_assays,
      flag_missing = frac_missing > max_missing
    ) |>
    arrange(desc(frac_missing), SampleID)

  if (any(assay_miss$n_missing < 0) || any(sample_miss$n_missing < 0)) {
    stop(sprintf("[%s] QC: negative missing counts -- duplicated rows detected. Data integrity failure.",
                 subcohort_name))
  }

  total_missing <- sum(assay_miss$n_missing)
  message(sprintf(
    "[%s] QC: %d assays x %d samples = %d cells; %d missing (%.2f%%); %d assays >%d%% missing; %d samples >%d%% missing",
    subcohort_name, n_assays, n_samples, n_assays * n_samples,
    total_missing, 100 * total_missing / (n_assays * n_samples),
    sum(assay_miss$flag_missing), round(max_missing * 100),
    sum(sample_miss$flag_missing), round(max_missing * 100)
  ))

  # --- 2. PCA-distance outlier screen ----------------------------------------
  pca_path <- file.path(out_dir, "pca.rds")
  if (!file.exists(pca_path)) {
    stop(sprintf("[%s] QC outlier screen needs %s -- run R/run_pca.R first.",
                 subcohort_name, pca_path))
  }
  scores <- readRDS(pca_path)$result

  pcs <- intersect(paste0("PC", seq_len(n_pcs)), names(scores))

  outliers <- scores |>
    select(SampleID, Individual, Group, all_of(pcs)) |>
    group_by(Group) |>
    group_modify(~ {
      m   <- as.matrix(.x[, pcs, drop = FALSE])
      ctr <- colMeans(m)
      d   <- sqrt(rowSums(sweep(m, 2, ctr, "-")^2))
      sc  <- mad(d)               # consistent sigma estimator
      if (sc <= 0) sc <- sd(d)    # tiny/degenerate group fallback
      z  <- if (sc > 0) (d - median(d)) / sc else rep(0, length(d))
      .x |> mutate(dist = d, z = z)
    }) |>
    ungroup() |>
    mutate(flagged = z > outlier_z) |>   # one-sided: distance can only be too large
    select(SampleID, Individual, Group, dist, z, flagged) |>
    arrange(desc(z))

  message(sprintf(
    "[%s] QC: outlier screen (distance to group centroid, PC1-%s, robust z > %g): %d of %d samples flagged",
    subcohort_name, tail(pcs, 1), outlier_z,
    sum(outliers$flagged), n_samples
  ))

  # --- 3. Plate x Group balance -----------------------------------------------
  plate_summary <- df |>
    distinct(SampleID, Group, Plate) |>
    count(Plate, Group, name = "n")

  n_plates <- n_distinct(plate_summary$Plate)
  biggest  <- plate_summary |>
    group_by(Plate) |>
    summarise(total = sum(n), .groups = "drop") |>
    slice_max(total, n = 1)

  message(sprintf(
    "[%s] QC: %d plates; largest %s (%d samples); plate x group counts in artifact",
    subcohort_name, n_plates, biggest$Plate, biggest$total
  ))

  # --- 4. Artifact ---------------------------------------------------------------
  params <- list(
    subcohort                 = subcohort_name,
    n_samples                 = n_samples,
    n_assays                  = n_assays,
    possible_cells            = n_assays * n_samples,
    n_missing_cells           = total_missing,
    frac_missing_cells        = total_missing / (n_assays * n_samples),
    missing_encoding          = "absent source rows; per-unit counts = expected - observed",
    max_missing_threshold     = max_missing,
    n_assays_flagged_missing  = sum(assay_miss$flag_missing),
    n_samples_flagged_missing = sum(sample_miss$flag_missing),
    outlier_screen            = list(
      method      = paste0("Euclidean distance to own-group centroid in PC1-",
                           tail(pcs, 1), ", robust z (median/MAD), one-sided"),
      z_threshold = outlier_z,
      n_flagged   = sum(outliers$flagged)
    ),
    n_plates                  = n_plates,
    lod_note                  = paste0("LOD flags not computable: LOD fields ",
                                       "absent from the harmonized dataset."),
    exclusion_policy          = paste0("QC flags informational; no samples ",
                                       "excluded (no pre-specified rule)."),
    finished_utc              = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  saveRDS(
    list(
      result = list(
        assay_missingness   = assay_miss,
        sample_missingness = sample_miss,
        outliers           = outliers,
        plate_summary      = plate_summary
      ),
      params = params
    ),
    file.path(out_dir, "qc.rds")
  )
  invisible(list(result = params, params = params))
}
