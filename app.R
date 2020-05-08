library(shiny)
library(plotly)
library(ggdendro)
library(igraph)
library(igraphdata)
library(gridExtra)
source("internal_funcs.R")
library(tictoc)
library(networkD3)
library(htmlwidgets)

ui <- fluidPage(
  titlePanel("AutoFocus Results"),

  selectInput(inputId = "outcome", label = "Select Disease Phenotype",
              choices = c("No outcome", colnames(R$samples)),
              selected = "No outcome"), #Outcome
  #sliderInput(inputId = "height", label = "Tree Height", min = 0, max = 1, value = 0.5),
  
      fluidRow(mainPanel("Module Browser", plotlyOutput("dendro"), plotlyOutput("clust"))),
      fluidRow( plotlyOutput("barplots")),
)

server <- function(input, output) {
  selected_node = reactiveValues(n = NA)
  observeEvent(event_data("plotly_click"),{
    d<-event_data("plotly_click")
    print(R$color[R$order_coords$index[which(R$order_coords$y ==d$y)]])
    if (R$color[R$order_coords$index[which(R$order_coords$y ==d$y)]] == "purple"){
      selected_node$n = R$order_coords$index[which(R$order_coords$y ==d$y)]
    }
  })
  
  output$dendro <- renderPlotly({
    
    y_axis <- list(
      title = "",
      showline = TRUE,
      showticklabels = F,
      mirror="ticks",
      showgrid = F,
      zeroline = F,
      linecolor = toRGB("black"),
      linewidth = 2
    )
    x_axis <- list(
      title = "",
      showline = TRUE,
      showticklabels = F,
      mirror="ticks",
      showgrid = F,
      zeroline = F,
      linecolor = toRGB("black"),
      linewidth = 2
    )

    #
    # ### Dendrogram ###
    #
    dend_network <- R$dend[[1]]
    dend_edges <- R$dend[[2]]
    x_axis$range <- R$views[as.double(selected_node$n)][[1]][[1]]
    y_axis$range <- R$views[as.double(selected_node$n)][[1]][[2]]

    dend_network<-plotly::layout(
      dend_network,
      shapes = dend_edges,
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = F
    )

    dend_network
  })
    
    ### Subnetwork ###
  
output$clust <- renderPlotly({
      
  axis <- list(
    title = "",
    showline = TRUE,
    showticklabels = F,
    mirror="ticks",
    showgrid = F,
    zeroline = F,
    linecolor = toRGB("black"),
    linewidth = 2
  )
  if (is.na(selected_node$n)){
    net = plotly_empty() %>% add_annotations(text = "Select a Significant Cluster to View", showarrow = F) %>% layout(xaxis = axis, yaxis = axis)
  }
  else{
  
    net <- R$msts[as.double(selected_node$n)][[1]][[1]]
    edge_shapes <- R$msts[as.double(selected_node$n)][[1]][[2]]
    net<-plotly::layout(
      net,
      shapes = edge_shapes,
      xaxis = axis,
      yaxis = axis,
      showlegend=F
    )
  }

  net
})
  
  
    
  ### Barplots ###
  output$barplots <- renderPlotly({
    axis <- list(
      
      showline = TRUE,
      showticklabels = F,
      mirror="ticks",
      showgrid = F,
      zeroline = F,
      linecolor = toRGB("black"),
      linewidth = 2
    )
    if (is.na(selected_node$n)){
      sunburst = plotly_empty() %>% 
        add_annotations(text = "Select a Significant Cluster for Annotation Details", showarrow = F) %>% 
        layout(xaxis = axis, yaxis = axis)

    }

    else{
      barplot_df <- data.frame(R$barplots[selected_node$n][[1]])
      sunburst <- plot_ly(barplot_df, labels = ~labels, parents = ~parents, values = ~values,type = 'sunburst', branchvalues = 'remainder')
    }
    
    sunburst %>% layout(
      title="Module Node Annotations")
    sunburst
    }
  )

    

}



shinyApp(ui, server)


