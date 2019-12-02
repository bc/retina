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
    radioButtons("side", "Eye", c("left","right")) ,
    tags$hr(),
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
      tags$hr()
    ),
    mainPanel(
      plotOutput("outline_indices_plot")
    )
  )
)

server <- function(input, output) {
  output$outline_indices_plot <- renderPlot({    

    roi_file <- input$outline_roi
    measurements_file <- input$datapoints

    if (!is.null(roi_file) & !is.null(measurements_file)){
    outline_indices_plot(load_roi(file.path(roi_file$datapath)), fread(file.path(measurements_file$datapath)))
    } else {
      plot(1:100)
    }

    # polygon(read.csv(falciform_file))
    
  },width = 1000, height = 1000)
}

shinyApp(ui, server)
