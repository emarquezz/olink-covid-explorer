# app.R -- Step 8: add Volcano tab (DE results explorer)
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)

npx <- readRDS("data/harmonized_npx.rds")
subcohort_choices <- sort(unique(npx$subcohort))

ui <- page_sidebar(
  title = "Olink COVID-19 Explorer",

  sidebar = sidebar(
    selectInput("subcohort", "Subcohort (sample matrix)",
                choices = subcohort_choices,
                selected = subcohort_choices[1]),
    selectInput("group",   "Sample group",     choices = NULL),
    selectInput("protein", "Protein (Assay)", choices = NULL)
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
    )
  )
)

server <- function(input, output, session) {

  # ---- 1. Cascading filters ------------------------------------------------
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
  # NOTE: plotly::validate masks shiny::validate (library load order), so
  # validation must be namespaced: shiny::validate(shiny::need(...)).
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
    d   <- volcano_data()
    pal <- c("Up in COVID"           = "#d73027",
             "Down in COVID"         = "#4575b4",
             "FDR-sig, small effect" = "#999999",
             "Not significant"       = "#d9d9d9")

    p <- ggplot(d, aes(x = logFC, y = -log10(adj.P.Val), color = category,
                       text = paste0(Assay, " (", UniProt, ")\n",
                                     "logFC = ", round(logFC, 2), "\n",
                                     "adj.P = ", format(adj.P.Val, digits = 2)))) +
      geom_point(alpha = 0.75, size = 2) +
      scale_color_manual(values = pal) +
      geom_vline(xintercept = c(-1, 1) * input$fc_thr,
                 linetype = "dashed", color = "grey45") +
      geom_hline(yintercept = -log10(input$fdr_thr),
                 linetype = "dashed", color = "grey45") +
      labs(x = "log2 fold change (POSITIVE - NEGATIVE)",
           y = "-log10(adj.P.Val)", color = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")

    ggplotly(p, tooltip = "text") |>
      layout(legend = list(orientation = "h"))
  })
}

shinyApp(ui, server)
