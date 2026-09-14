# app.R -- Step 6: subcohort dropdown + cascading filters + summary cards
library(shiny)
library(bslib)
library(DT)
library(dplyr)

npx <- readRDS("data/harmonized_npx.rds")

subcohort_choices <- sort(unique(npx$subcohort))

ui <- page_sidebar(
  title = "Olink COVID-19 Explorer",

  sidebar = sidebar(
    selectInput("subcohort", "Subcohort (sample matrix)",
                choices = subcohort_choices,
                selected = subcohort_choices[1]),
    selectInput("group",     "Sample group",      choices = NULL),
    selectInput("protein",   "Protein (Assay)",   choices = NULL)
  ),

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
)

server <- function(input, output, session) {

  # ---- 1. Cascading filters ----------------------------------------------
  # Fires on app startup too (observeEvent's ignoreInit defaults to FALSE,
  # and input$subcohort already has a value) -- that's what populates the
  # two NULL-choices dropdowns on first load.
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

  # ---- 2. Rows for the current subcohort (+ group if not "All") ----------
  current <- reactive({
    req(input$subcohort, input$group)          # don't run until both are set
    d <- npx |> filter(subcohort == input$subcohort)
    if (input$group != "All") d <- d |> filter(Group == input$group)
    d
  })

  # ---- 3. Rows for the selected protein -----------------------------------
  filtered <- reactive({
    req(input$protein)
    current() |> filter(Assay == input$protein)
  })

  # ---- 4. Outputs ---------------------------------------------------------
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
}

shinyApp(ui, server)
