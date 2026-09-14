# R/load_data.R
# Harmonization layer for the Gisby et al. (2021) Olink COVID-19 dataset.
#
# Verified schema (2026-09-14, from the raw CSVs):
#   NPX files (long, 8 cols): SampleID, Individual_ID, UniProt, GeneID,
#     Assay, NPX, Panel, Index   (Index = run-sheet position; dropped)
#   Plasma sample-level: 23 cols (severity, labs, time columns, ...)
#   Serum sample-level:  7 cols (NO severity/labs/time) -> optional mapping
#   Group levels verified: POSITIVE / NEGATIVE
#   Age is provided pre-binned (e.g. "(60,80]") and kept as character.
#   No LOD columns exist in this public release; missingness is encoded as
#   ABSENT ROWS (not NA) and is handled downstream (QC/DE steps).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

# --- Canonical schema -------------------------------------------------------
# Left = canonical name used everywhere downstream; right = raw CSV column.
canonical_npx_cols <- c(
  SampleID   = "SampleID",
  Individual = "Individual_ID",
  UniProt    = "UniProt",     # cross-platform harmonization key
  GeneID     = "GeneID",
  Assay      = "Assay",
  Panel      = "Panel",
  NPX        = "NPX"
)

canonical_meta_required <- c(
  SampleID   = "SampleID",
  Individual = "Individual_ID",
  Plate      = "Plate_ID",     # kept: enables plate/batch QC downstream
  Group      = "Case_Control",
  Sex        = "Sex",
  Age        = "Age"           # pre-binned in source; kept as character
)

# Present in plasma, ABSENT in serum -> filled with NA, never a hard error.
canonical_meta_optional <- c(
  Severity = "WHO_Severity_Peak"
)

EXPECTED_GROUPS <- c("POSITIVE", "NEGATIVE")

# Rename per mapping. `optional` = raw column names that may be absent
# (filled with NA + logged) instead of stopping the pipeline.
# All other mapped columns are required; any absence is a hard error.
standardize <- function(df, mapping, label, optional = character(0)) {
  raw_cols    <- unname(mapping)
  missing_req <- setdiff(setdiff(raw_cols, optional), names(df))
  if (length(missing_req) > 0) {
    stop(sprintf(
      "[%s] Missing expected columns: %s\n  Available: %s",
      label, paste(missing_req, collapse = ", "),
      paste(names(df), collapse = ", ")
    ))
  }

  missing_opt <- setdiff(optional, names(df))
  if (length(missing_opt) > 0) {
    for (col in missing_opt) df[[col]] <- NA_character_
    message(sprintf("[%s] Optional column(s) absent, filled with NA: %s",
                    label, paste(missing_opt, collapse = ", ")))
  }

  df |>
    select(all_of(raw_cols)) |>
    rename(!!!setNames(raw_cols, names(mapping)))
}

load_subcohort <- function(npx_path, meta_path, subcohort) {
  npx  <- read_csv(npx_path,  show_col_types = FALSE)
  meta <- read_csv(meta_path, show_col_types = FALSE)

  npx <- standardize(npx, canonical_npx_cols, paste0(subcohort, "/npx"))

  # One call for the whole metadata schema: required columns must exist,
  # optional ones (severity, absent in serum) are NA-filled and logged.
  meta <- standardize(meta,
                      c(canonical_meta_required, canonical_meta_optional),
                      paste0(subcohort, "/meta"),
                      optional = unname(canonical_meta_optional))

  # Verified-constants assertion: unexpected group levels stop the pipeline
  # rather than silently flowing into the design matrix downstream.
  bad <- setdiff(unique(meta$Group), EXPECTED_GROUPS)
  if (length(bad) > 0) {
    stop(sprintf("[%s/meta] Unexpected Case_Control levels: %s",
                 subcohort, paste(bad, collapse = ", ")))
  }

  # Join on BOTH keys. This doubles as a referential-integrity check:
  # a SampleID present in both files but with conflicting Individual_ID
  # fails to match and is counted as a key mismatch below.
  n_before <- nrow(npx)
  df <- npx |>
    left_join(meta, by = c("SampleID", "Individual")) |>
    mutate(subcohort = subcohort)

  if (nrow(df) != n_before) {
    stop(sprintf(
      "[%s] Join changed row count (%d -> %d): duplicated keys in metadata.",
      subcohort, n_before, nrow(df)
    ))
  }

  orphans      <- length(setdiff(npx$SampleID, meta$SampleID))
  key_mismatch <- df |>
    filter(SampleID %in% meta$SampleID, is.na(Group)) |>
    nrow()
  unused       <- length(setdiff(meta$SampleID, npx$SampleID))

  n_na <- sum(is.na(df$Group))
  if (n_na > 0) {
    df <- filter(df, !is.na(Group))
  }
  message(sprintf(
    paste0("[%s] %d NPX rows; dropped %d without Group (%d orphan samples, ",
           "%d key mismatches); %d meta samples have no NPX data"),
    subcohort, nrow(df), n_na, orphans, key_mismatch, unused
  ))

  df
}

# Main entry point: run once to produce data/harmonized_npx.rds
harmonize <- function(data_dir = here("data"),
                      out_path = NULL) {
  if (is.null(out_path)) out_path <- file.path(data_dir, "harmonized_npx.rds")

  plasma <- load_subcohort(
    file.path(data_dir, "plasma_npx_level.csv"),
    file.path(data_dir, "plasma_sample_level.csv"),
    "plasma"
  )
  serum <- load_subcohort(
    file.path(data_dir, "serum_npx_level.csv"),
    file.path(data_dir, "serum_sample_level.csv"),
    "serum"
  )

  combined <- bind_rows(plasma, serum)

  # Cross-matrix overlap: the cross-study story in numbers.
  p_up <- unique(plasma$UniProt)
  s_up <- unique(serum$UniProt)
  message(sprintf(
    paste0("Harmonized: %d plasma rows, %d serum rows | ",
           "%d proteins shared (plasma %d, serum %d) | ",
           "%d shared SampleIDs | %d shared individuals"),
    nrow(plasma), nrow(serum),
    length(intersect(p_up, s_up)), length(p_up), length(s_up),
    length(intersect(unique(plasma$SampleID), unique(serum$SampleID))),
    length(intersect(unique(plasma$Individual), unique(serum$Individual)))
  ))

  saveRDS(combined, out_path)
  invisible(combined)
}
