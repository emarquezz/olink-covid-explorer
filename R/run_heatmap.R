# R/run_heatmap.R
# Top-DE protein heatmaps per subcohort (ComplexHeatmap -> PNG artifacts).
#
# Renders TWO figures per subcohort:
#   1. heatmap.png             -- per-sample view: top-N proteins x all
#      samples, columns ordered by Group -> Severity (clinical order) ->
#      SampleID, NOT clustered (blocks stay contiguous). Group/Severity
#      column annotations + DE-direction row bar. The QC/structure view.
#   2. heatmap_categories.png  -- category-average view: top-N proteins x
#      Group/Severity categories (plasma: NEGATIVE, mild, moderate, severe,
#      critical; serum falls back to Group: NEGATIVE, POSITIVE), values =
#      mean per-protein z-score within category, columns in clinical order
#      (dose-response reads left -> right), NOT clustered. The story view.
#
# Design notes:
#   - Same matrix + missingness recipe as run_de / run_pca.
#   - Row selection + DE-direction bar come from
#     outputs/<subcohort>/de_results.rds (pipeline consuming another
#     pipeline's artifact).
#   - Row-wise z-scores: rows ARE proteins here, so t(scale(t(mat))) is the
#     CORRECT form (in run_pca it was scale(t(mat)) because samples had to
#     be rows). Same operation family, different target.
#   - Category figure clusters rows on the category-AVERAGE matrix, so
#     proteins group by severity profile; its row order may differ from
#     the per-sample figure. Documented in params, intentional.
#   - Per-category PER-SAMPLE heatmaps deliberately NOT rendered: severity
#     bins have small n, within-bin variance is noise-dominated once the
#     COVID signal is conditioned away, and cross-figure comparison is
#     poor. The per-sample figure already carries that structure; the
#     category figure is its compact summary.
#
# Reads:  data/harmonized_npx.rds (via caller),
#         outputs/<subcohort>/de_results.rds
# Writes: outputs/<subcohort>/heatmap.png,
#         outputs/<subcohort>/heatmap_categories.png,
#         outputs/<subcohort>/heatmap.rds (top table + params)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(here)
})

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  stop("ComplexHeatmap not installed. Run: renv::install('bioc::ComplexHeatmap')")
}

run_heatmap <- function(npx, subcohort_name, out_dir,
                        top_n = 50, max_missing = 0.20) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df <- npx |> filter(subcohort == .env$subcohort_name)

  # --- 0. DE artifact (pipeline prerequisite: hard error, not validate) ----
  de_path <- file.path(out_dir, "de_results.rds")
  if (!file.exists(de_path)) {
    stop(sprintf("[%s] heatmap needs %s -- run R/run_de.R first.",
                 subcohort_name, de_path))
  }
  de <- readRDS(de_path)

  # --- 1. Protein x sample matrix (same recipe as run_de/run_pca) ----------
  mat <- df |>
    select(Assay, SampleID, NPX) |>
    pivot_wider(names_from = SampleID, values_from = NPX) |>
    column_to_rownames("Assay") |>
    as.matrix()

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
    "[%s] heatmap matrix: %d proteins x %d samples (%d dropped >%d%% missing; %d cells median-imputed)",
    subcohort_name, nrow(mat), ncol(mat), n_dropped,
    round(max_missing * 100), n_imputed
  ))

  # --- 2. Top-N selection from the DE artifact -------------------------------
  top_tab <- de$result |>
    arrange(adj.P.Val, desc(abs(logFC))) |>
    slice_head(n = top_n) |>
    mutate(
      direction = factor(
        if_else(adj.P.Val < 0.05,
                if_else(logFC > 0, "Up in COVID", "Down in COVID"),
                "Not significant"),
        levels = c("Up in COVID", "Down in COVID", "Not significant")
      )
    )

  # --- 3. Row-wise z-scores, subset to top proteins --------------------------
  zmat <- t(scale(t(mat)))                       # rows ARE proteins here
  zmat <- zmat[rownames(zmat) %in% top_tab$Assay, , drop = FALSE]
  top_tab <- top_tab[match(rownames(zmat), top_tab$Assay), ]

  # --- 4. Column order: Group, then Severity, then SampleID -----------------
  clinical <- c("NEGATIVE", "mild", "moderate", "severe", "critical")
  ann <- df |>
    distinct(SampleID, Group, Severity) |>
    filter(SampleID %in% colnames(zmat)) |>
    mutate(
      Group    = factor(Group, levels = c("NEGATIVE", "POSITIVE")),
      sev_rank = match(Severity, clinical)
    ) |>
    arrange(Group, sev_rank, SampleID)

  ord  <- ann$SampleID
  zmat <- zmat[, ord, drop = FALSE]   # ann rows now align with zmat cols

  # --- 5. Annotations (per-sample figure) ------------------------------------
  group_col <- c(POSITIVE = "#d73027", NEGATIVE = "#4575b4")

  sev_levels <- intersect(clinical, unique(na.omit(ann$Severity)))
  show_sev   <- length(sev_levels) > 0
  ann_cols   <- list(Group = group_col)
  if (show_sev) {
    ramp_n  <- sum(sev_levels != "NEGATIVE")
    sev_col <- setNames(
      c(if ("NEGATIVE" %in% sev_levels) "#999999",
        if (ramp_n > 0)
          colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(ramp_n)),
      sev_levels
    )
    ann_cols$Severity <- sev_col
  }

  ann_df <- ann[, names(ann_cols), drop = FALSE]   # Group (+ Severity)
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = ann_df, col = ann_cols, na_col = "grey90"
  )

  dir_col <- c("Up in COVID"     = "#d73027",
               "Down in COVID"   = "#4575b4",
               "Not significant" = "#d9d9d9")
  left_annotation <- ComplexHeatmap::rowAnnotation(
    DE = top_tab$direction, col = list(DE = dir_col)
  )

  # --- 6. Per-sample heatmap -> heatmap.png ----------------------------------
  hm_col <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "RdBu")))(101)  # blue (low) -> red (high)

  ht <- ComplexHeatmap::Heatmap(
    zmat,
    name              = "NPX (z)",
    col               = hm_col,
    top_annotation    = top_annotation,
    left_annotation   = left_annotation,
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_column_names = FALSE,
    row_names_gp      = grid::gpar(fontsize = 7),
    column_title      = sprintf("%s: top %d proteins by DE evidence",
                                subcohort_name, nrow(zmat)),
    column_title_gp   = grid::gpar(fontsize = 12)
  )

  png_path <- file.path(out_dir, "heatmap.png")
  w_px <- max(1000, ncol(zmat) * 8 + 350)
  h_px <- max(700,  nrow(zmat) * 18 + 320)
  png(png_path, width = w_px, height = h_px)
  ComplexHeatmap::draw(ht)
  dev.off()

  # --- 7. Category-average heatmap -> heatmap_categories.png ----------------
  # Category = Severity when it carries values (plasma; controls carry the
  # literal "NEGATIVE"), else Group (serum: Severity all-NA by design).
  use_sev <- sum(!is.na(ann$Severity)) > 0
  ann$Category <- if (use_sev) ann$Severity else as.character(ann$Group)

  cat_levels <- if (use_sev) clinical else c("NEGATIVE", "POSITIVE")
  cats  <- intersect(cat_levels, unique(na.omit(ann$Category)))
  n_excl <- nrow(ann) - sum(ann$Category %in% cats)
  if (n_excl > 0) {
    message(sprintf(
      "[%s] category heatmap: %d samples excluded (NA/unknown category)",
      subcohort_name, n_excl
    ))
  }

  # Mean per-protein z-score within each category (compact summary of the
  # same data -- no new modeling). ann rows align with zmat columns.
  avg <- sapply(cats, function(cl) {
    rowMeans(zmat[, which(ann$Category == cl), drop = FALSE])
  })
  colnames(avg) <- cats

  if (use_sev) {
    ramp_n  <- sum(cats != "NEGATIVE")
    cat_pal <- setNames(
      c(if ("NEGATIVE" %in% cats) "#999999",
        if (ramp_n > 0)
          colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(ramp_n)),
      cats
    )
  } else {
    cat_pal <- c(NEGATIVE = "#4575b4", POSITIVE = "#d73027")
  }

  cat_ann <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(Category = factor(colnames(avg), levels = cats)),
    col = list(Category = cat_pal)
  )

  ht_cat <- ComplexHeatmap::Heatmap(
    avg,
    name              = "mean NPX (z)",
    col               = hm_col,
    top_annotation    = cat_ann,
    left_annotation   = left_annotation,
    cluster_columns   = FALSE,     # clinical order IS the story
    cluster_rows      = TRUE,      # on the category-average matrix
    show_column_names = TRUE,
    column_names_gp   = grid::gpar(fontsize = 10),
    row_names_gp      = grid::gpar(fontsize = 7),
    column_title      = sprintf("%s: top %d proteins, mean z per category",
                                subcohort_name, nrow(avg)),
    column_title_gp   = grid::gpar(fontsize = 12)
  )

  png2_path <- file.path(out_dir, "heatmap_categories.png")
  w2_px <- 700
  h2_px <- max(700, nrow(zmat) * 18 + 320)
  png(png2_path, width = w2_px, height = h2_px)
  ComplexHeatmap::draw(ht_cat)
  dev.off()

  counts_str <- paste(
    sprintf("%s=%d", cats,
            as.integer(table(factor(ann$Category, levels = cats)))),
    collapse = ", "
  )
  message(sprintf(
    "[%s] category heatmap: %s (%d excluded) -> heatmap_categories.png",
    subcohort_name, counts_str, n_excl
  ))

  # --- 8. Artifact -------------------------------------------------------------
  params <- list(
    subcohort          = subcohort_name,
    top_n              = nrow(zmat),
    selection_rule     = "adj.P.Val ascending, tie-break |logFC| descending",
    n_up               = sum(top_tab$direction == "Up in COVID"),
    n_down             = sum(top_tab$direction == "Down in COVID"),
    n_samples          = ncol(zmat),
    column_order       = paste0("Group (NEGATIVE first), Severity (clinical ",
                                "order), SampleID; not clustered"),
    row_clustering     = "Euclidean; rows only",
    annotations        = names(ann_cols),
    max_missing        = max_missing,
    n_proteins_dropped = n_dropped,
    n_cells_imputed    = n_imputed,
    png                = "heatmap.png",
    png_dims           = c(width = w_px, height = h_px),
    category_view      = list(
      categories      = cats,
      n_per_category  = table(factor(ann$Category, levels = cats)),
      n_excluded      = n_excl,
      value           = "mean per-protein z-score within category",
      category_source = if (use_sev) "Severity" else "Group",
      row_clustering  = paste0("Euclidean on category-average matrix; row ",
                               "order may differ from per-sample figure"),
      png             = "heatmap_categories.png",
      png_dims        = c(width = w2_px, height = h2_px)
    ),
    finished_utc       = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  saveRDS(list(result = top_tab, params = params),
          file.path(out_dir, "heatmap.rds"))
  message(sprintf(
    "[%s] heatmap: top %d by adj.P (%d up / %d down) x %d samples -> heatmap.png + heatmap_categories.png",
    subcohort_name, nrow(zmat), params$n_up, params$n_down, ncol(zmat)
  ))
  invisible(list(result = top_tab, params = params))
}
