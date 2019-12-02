source("helper_functions.r")
library(shiny)
# roi_path <- "sample_retina/outline.roi"

# measurements <- fread("sample_retina/datapoints.csv")
# values <- fread("sample_retina/xyz.csv")[,3] %>% unlist %>% as.numeric
# measurements[,measurement := values]

# falciform_coords_raw <- fread("sample_retina/falciform.csv")

# save_outline_indices_plot(roi_path,measurements,"output/outline_coordinates.pdf")

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
    radioButtons("side", "eye", c("left","right")) ,
      fileInput("datapoints", "Choose CSV File for Datapoints",
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv")
        ),
      fileInput("outline_roi", "Choose outline ROI File",
        accept = c(
          "text/roi",
          ".roi")
        ),
       fileInput("falciform", "Choose falciform/optic disc CSV File",
        accept = c(
          "text/roi",
          ".roi")
        ),
      tags$hr(),
      checkboxInput("header", "Header", TRUE)
    ),
    mainPanel(
      tableOutput("table_contents")
    )
  )
)

server <- function(input, output) {
  output$table_contents <- renderTable({
    # input$datapoints will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$datapoints

    if (is.null(inFile))
      return(NULL)

    read.csv(inFile$datapath, header = input$header)
  })
  output$outline_indices_plot <- renderPlot({
    roi_file <- input$outline_roi

    if (is.null(roi_file))
      return(NULL)

    measurements_file <- input$datapoints

    if (is.null(measurements_file))
      return(NULL)

    falciform_file <- input$falciform
    if (is.null(falciform_file))
      return(NULL)

    outline_indices_plot(roi_file, measurements_file)
    polygon(read.csv(falciform_file))
    
  })
}

shinyApp(ui, server)
