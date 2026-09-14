# app.R -- Step 5: minimal app: one protein dropdown, one table.
library(shiny)
library(bslib)
library(DT)
library(dplyr)

# Load harmonized data ONCE, at top level (outside ui/server).
# Shiny sets the working directory to the app folder, so the
# relative path resolves regardless of where you launch from.
npx <- readRDS("data/harmonized_npx.rds")

assays <- sort(unique(npx$Assay))
# Default to a biologically interesting protein if present.
default_assay <- if ("IL6" %in% assays) "IL6" else assays[1]

# ---------------------------------------------------------------------------
# UI: what the page looks like
# ---------------------------------------------------------------------------
ui <- page_sidebar(
  title = "Olink COVID-19 Explorer",

  sidebar = sidebar(
    selectInput(
      inputId  = "protein",
      label    = "Protein (Assay)",
      choices  = assays,          # selectInput is searchable: just type
      selected = default_assay
    )
  ),

  card(
    card_header("NPX values for selected protein"),
    DTOutput("npx_table")
  )
)

# ---------------------------------------------------------------------------
# Server: reactive logic
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # A reactive expression: caches its result and reruns ONLY when
  # input$protein changes. Assay/UniProt/Panel are constant for one
  # protein, so we drop them from the display and show what varies.
  filtered <- reactive({
    npx |>
      filter(Assay == input$protein) |>
      select(SampleID, Individual, NPX, Group, Sex, Age,
             Severity, Plate, subcohort)
  })

  # Fill the DTOutput slot. Calling filtered() (with parentheses)
  # is what creates the dependency edge.
  output$npx_table <- renderDT({
    datatable(
      filtered(),
      options  = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)
