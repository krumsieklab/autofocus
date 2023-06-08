#' Launch shiny app
#'
#' Main autofocus function to launch shiny app.
#'
#' @param R R struct
#' @param anno_list Vector of annotations to extract from R struct
#'
#' @return None (launches app)
#'
#' @import shiny
#' @import shinydashboard
#' @import plotly
#' @import networkD3
#'
#' @export
run_autofocus <- function(R, anno_list = c("")){

  # move the anno_list argument to the intialize functino

  appDir <<- system.file("shiny-app", "autofocus", package = "autofocus")

  if(appDir == "") stop("Shiny app directory not found. Try re-installing autofocus package.")

  R$annos <<- R$annos[anno_list]
  color_palette<<- MetBrewer::met.brewer("Hiroshige",10)
  body <- shinydashboard::dashboardBody(
    fluidRow(tabBox(title=NULL,
                    width=12,
                    tabPanel("Tree View",plotly::plotlyOutput("dendro")),
                    tabPanel("Peak List",DT::dataTableOutput("all_modules_table")),
                    tabPanel("Analyte List",DT::dataTableOutput("analyte_table")),
                    tabPanel("Peak Annotations",reactable::reactableOutput("anno_table")),
                    sliderInput("threshold", "Density Threshold:", min=0, max=1, value=0.5)),
             tabBox(title = NULL,
                    width=12,
                    tabPanel("Module members", DT::dataTableOutput("single_module_table")),
                    tabPanel("Module Network and Drivers",
                             fluidRow(
                               column(networkD3::forceNetworkOutput("network"),width=9),
                               column(DT::dataTableOutput("drivers"),width=11)))))
  )

  #source(paste0(appDir, '/ui.R'))

  ui<-dashboardPage(dashboardHeader(title="AutoFocus Run Results"),
                    dashboardSidebar(disable = TRUE),
                    body)

  source(paste0(appDir, '/server.R'))

  shiny::runApp(list(ui = ui, server = server),
                launch.browser = TRUE)

}
