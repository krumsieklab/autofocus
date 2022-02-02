library(plotly)
library(shiny)

ui <- fluidPage(
  plotlyOutput("scatter")
)
server <- function(input, output) {
  
  output$scatter <- renderPlotly({
    data <- data.frame(x=sample.int(10000,10000), y = sample.int(10000,10000))
    plot_ly(data, x = ~x, y = ~y) %>% layout(xaxis=x_axis, yaxis=y_axis)
  })
  
  scatterProxy <- plotlyProxy("scatter")
  
  observeEvent(event_data("plotly_click"),{
    d<-event_data("plotly_click")
    xrange <- c((d$x- 1),(d$x+ 1))
    yrange <- c((d$y- 1),(d$y- 1))
    plotlyProxyInvoke(scatterProxy, "relayout", list(xaxis = list(range = xrange), yaxis = list(range = yrange)))
  })
  

}
shinyApp(ui, server)