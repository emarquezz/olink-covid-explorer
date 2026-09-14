# app.R -- Steps 5-10: data table + volcano + PCA tabs
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(plotly)

npx <- readRDS("data/harmonized_npx.rds")
subcohort_choices <- sort(unique(npx$subcohort))

ui <- page_sidebar(
  title = "Olink COVID-19 Explorer",

  sidebar = sidebar(
    selectInput("subcohort", "Subcohort (sample matrix)",
                choices = subcohort_choices,
                selected = subcohort_choices[1]),
    selectInput("group",   "Sample group (Data table & PCA)", choices = NULL),
    selectInput("protein", "Protein (Assay)",                 choices = NULL)
  ),

  navset_card_tab(
    full_screen = TRUE,
    id = "tabs",

    # ---- Tab 1: Data table ------------------------------------------------
    nav_panel("Data table", value = "table",
              layout_columns(
                col_widths = c(6, 6),
                card(card_header("Samples in current selection"),
                     textOutput("n_samples")),
                card(card_header("Proteins in current selection"),
                     textOutput("n_proteins"))
              ),
              card(
                card_header("NPX values for selected protein"),
                DTOutput("npx_table")
              )
    ),

    # ---- Tab 2: Volcano ----------------------------------------------------
    nav_panel("Volcano", value = "volcano",
              card(
                card_header("Differential expression: COVID POSITIVE vs NEGATIVE"),
                p("Cohort-level statistic computed by R/run_de.R on the full subcohort; sidebar sample filters do not apply to this tab.",
                  class = "text-muted", style = "font-size: 0.85em; margin: 0 1rem 0.5rem;"),
                layout_columns(
                  col_widths = c(6, 6),
                  sliderInput("fc_thr", "|logFC| threshold",
                              min = 0, max = 2, value = 0.5, step = 0.1),
                  sliderInput("fdr_thr", "FDR threshold (adj.P)",
                              min = 0.001, max = 0.25, value = 0.05, step = 0.005)
                ),
                textOutput("volcano_counts"),
                plotlyOutput("volcano", height = "500px")
              )
    ),

    # ---- Tab 3: PCA --------------------------------------------------------
    nav_panel("PCA", value = "pca",
              card(
                card_header("Sample PCA: PC1 vs PC2"),
                selectInput("pca_color", "Color by", choices = "Group"),
                textOutput("pca_subtitle"),
                plotlyOutput("pca_plot", height = "520px")
              )
    )
  )
)

server <- function(input, output, session) {

  # ---- 1. Cascading filters (Data table + PCA point filtering) ------------
  observeEvent(input$subcohort, {
    df <- npx |> filter(subcohort == input$subcohort)

    groups <- sort(unique(df$Group))
    updateSelectInput(session, "group",
                      choices = c("All", groups),
                      selected = "All")

    proteins <- sort(unique(df$Assay))
    default  <- if ("IL6" %in% proteins) "IL6" else proteins[1]
    updateSelectInput(session, "protein",
                      choices = proteins,
                      selected = default)
  })

  current <- reactive({
    req(input$subcohort, input$group)
    d <- npx |> filter(subcohort == input$subcohort)
    if (input$group != "All") d <- d |> filter(Group == input$group)
    d
  })

  filtered <- reactive({
    req(input$protein)
    current() |> filter(Assay == input$protein)
  })

  output$n_samples  <- renderText(n_distinct(current()$SampleID))
  output$n_proteins <- renderText(n_distinct(current()$Assay))

  output$npx_table <- renderDT({
    datatable(
      filtered() |>
        select(SampleID, Individual, NPX, Group, Sex, Age,
               Severity, Plate, subcohort),
      options  = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # ---- 2. Volcano -----------------------------------------------------------
  # plotly::validate masks shiny::validate (library load order), so
  # validation calls are namespaced: shiny::validate(shiny::need(...)).
  de_results <- reactive({
    req(input$subcohort)
    path <- file.path("outputs", input$subcohort, "de_results.rds")
    shiny::validate(shiny::need(
      file.exists(path),
      "DE results not found -- run R/run_de.R first."
    ))
    readRDS(path)$result
  })

  volcano_data <- reactive({
    de_results() |>
      mutate(
        category = case_when(
          adj.P.Val < input$fdr_thr & logFC >  input$fc_thr ~ "Up in COVID",
          adj.P.Val < input$fdr_thr & logFC < -input$fc_thr ~ "Down in COVID",
          adj.P.Val < input$fdr_thr                          ~ "FDR-sig, small effect",
          TRUE                                                ~ "Not significant"
        )
      )
  })

  output$volcano_counts <- renderText({
    d <- volcano_data()
    sprintf("%d up / %d down in COVID (|logFC| > %.1f, adj.P < %.3f) of %d proteins",
            sum(d$category == "Up in COVID"),
            sum(d$category == "Down in COVID"),
            input$fc_thr, input$fdr_thr, nrow(d))
  })

  output$volcano <- renderPlotly({
    pal <- c("Up in COVID"           = "#d73027",
             "Down in COVID"         = "#4575b4",
             "FDR-sig, small effect" = "#999999",
             "Not significant"       = "#d9d9d9")

    d <- volcano_data() |>
      mutate(
        category = factor(category, levels = names(pal)),
        p_neg    = -log10(pmax(adj.P.Val, 1e-300))   # Inf guard
      )

    ymax <- max(d$p_neg) * 1.05
    xmin <- min(d$logFC) * 1.05
    xmax <- max(d$logFC) * 1.05

    # Threshold lines: 2 vertical + 1 horizontal, one trace
    thr <- data.frame(
      x    = c( input$fc_thr, -input$fc_thr, xmin),
      xend = c( input$fc_thr, -input$fc_thr, xmax),
      y    = c(0, 0, -log10(input$fdr_thr)),
      yend = c(ymax, ymax, -log10(input$fdr_thr))
    )

    plot_ly(d, x = ~logFC, y = ~p_neg,
            color = ~category, colors = pal,
            type = "scatter", mode = "markers",
            marker = list(size = 7, opacity = 0.75),
            text = ~paste0("<b>", Assay, "</b> (", UniProt, ")<br>",
                           "logFC: ", round(logFC, 2),
                           " | adj.P: ", format(adj.P.Val, digits = 2)),
            hoverinfo = "text") |>
      add_segments(data = thr, x = ~x, xend = ~xend, y = ~y, yend = ~yend,
                   inherit = FALSE, showlegend = FALSE,
                   line = list(color = "grey45", dash = "dash")) |>
      layout(
        xaxis = list(title = "log2 fold change (POSITIVE - NEGATIVE)",
                     zeroline = FALSE),
        yaxis = list(title = "-log10(adj.P.Val)", zeroline = FALSE),
        legend = list(orientation = "h", x = 0, y = 1.08),
        margin = list(t = 80, r = 20)
      )
  })

  # ---- 3. PCA ---------------------------------------------------------------
  pca_data <- reactive({
    req(input$subcohort)
    path <- file.path("outputs", input$subcohort, "pca.rds")
    shiny::validate(shiny::need(
      file.exists(path),
      "PCA results not found -- run R/run_pca.R first."
    ))
    readRDS(path)
  })

  # Adaptive color-by: choices computed FROM THE DATA, not hardcoded per
  # subcohort. A column qualifies only if it exists and has >1 distinct
  # value -- so serum (Severity all-NA) and single-valued columns drop out.
  observeEvent(input$subcohort, {
    res <- pca_data()$result
    candidates <- c("Group", "Sex", "Age", "Severity", "Plate")
    ok <- candidates[vapply(candidates, function(v) {
      v %in% names(res) && n_distinct(res[[v]], na.rm = TRUE) > 1
    }, logical(1))]
    updateSelectInput(session, "pca_color",
                      choices = ok, selected = "Group")
  })

  # PCA points respect the sidebar group filter. Statistically safe:
  # PC axes come from the full-subcohort artifact; filtering only
  # selects which samples are displayed on those fixed axes.
  pca_filtered <- reactive({
    req(input$group)
    d <- pca_data()$result
    if (input$group != "All") d <- d |> filter(Group == input$group)
    d
  })

  output$pca_subtitle <- renderText({
    p <- pca_data()
    sprintf("PC1 %.1f%% | PC2 %.1f%% | %d of %d samples | colored by %s",
            p$params$var_explained[1] * 100,
            p$params$var_explained[2] * 100,
            nrow(pca_filtered()), nrow(p$result),
            input$pca_color)
  })

  output$pca_plot <- renderPlotly({
    req(input$pca_color)
    d   <- pca_filtered()
    col <- input$pca_color

    d <- d |>
      mutate(
        color_by = factor(.data[[col]]),
        tip = paste0("<b>", SampleID, "</b> (", Individual, ")<br>",
                     col, ": ", .data[[col]], "<br>",
                     "PC1: ", round(PC1, 1), " | PC2: ", round(PC2, 1))
      )

    # --- Color strategy: fixed pair for Group; ORDINAL gradients for
    # Severity / Age (darker = more severe / older); qualitative otherwise.
    if (col == "Group") {
      d <- d |> mutate(color_by = factor(color_by,
                                         levels = c("POSITIVE", "NEGATIVE")))
      pal <- c(POSITIVE = "#d73027", NEGATIVE = "#4575b4")

    } else if (col == "Severity") {
      clinical <- c("NEGATIVE", "mild", "moderate", "severe", "critical")
      present  <- intersect(clinical, levels(droplevels(d$color_by)))
      ramp_n   <- length(present) - as.integer("NEGATIVE" %in% present)
      ramp     <- colorRampPalette(
        RColorBrewer::brewer.pal(9, "YlOrRd"))(ramp_n)
      pal <- setNames(
        c(if ("NEGATIVE" %in% present) "#999999", ramp),
        present
      )
      d <- d |> mutate(color_by = factor(color_by, levels = present))

    } else if (col == "Age") {
      # order bins by numeric lower bound: (0,20] < (20,40] < ...
      # ramp anchored at brewer stops 3-9 so the youngest bin stays visible
      lv   <- levels(d$color_by)
      lv   <- lv[order(as.numeric(sub("^\\((\\d+).*$", "\\1", lv)))]
      ramp <- colorRampPalette(
        RColorBrewer::brewer.pal(9, "Blues")[5:9])(length(lv))
      pal  <- setNames(ramp, lv)
      d    <- d |> mutate(color_by = factor(color_by, levels = lv))

    } else {
      pal <- setNames(scales::hue_pal()(nlevels(d$color_by)),
                      levels(d$color_by))
    }

    plot_ly(d, x = ~PC1, y = ~PC2,
            color = ~color_by, colors = pal,
            type = "scatter", mode = "markers",
            marker = list(size = 8, opacity = 0.8,
                          line = list(width = 0.5, color = "white")),
            text = ~tip, hoverinfo = "text") |>
      layout(
        xaxis = list(title = "PC1", zeroline = FALSE),
        yaxis = list(title = "PC2", zeroline = FALSE),
        legend = list(title = list(text = col),
                      orientation = "h", x = 0, y = 1.08),
        margin = list(t = 80, r = 20)
      )
  })
}

shinyApp(ui, server)
