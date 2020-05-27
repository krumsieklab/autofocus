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
library(dendextend)



ui <- fluidPage(
  titlePanel("AutoFocus Results"),

  selectInput(inputId = "outcome", label = "Select Disease Phenotype",
              choices = c("No outcome", colnames(R$samples)),
              selected = "No outcome"), #Outcome

      fluidRow(mainPanel("Module Browser", plotlyOutput("dendro"), plotlyOutput("clust"))),
      fluidRow( plotlyOutput("barplots")),DT::dataTableOutput("mytable")
)


server <- function(input, output) {
  dend_xy <- R$HCL %>% as.dendrogram %>%get_nodes_xy()
  order_coords <- data.frame(index = R$order, 
                             x = round(dend_xy[,1], digits = 10), 
                             y = round(dend_xy[,2], digits = 10)) 
  
  order_coords <- order_coords[order(order_coords$index),]
  
  
  selected_node = reactiveValues(n = NA)
  observeEvent(event_data("plotly_click"),{
    d<-event_data("plotly_click")
    if (d$y != 0){
      selected_node$n = order_coords$index[which(order_coords$y ==d$y)]
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
    
    view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order) 
    x <- c((min(view_range)-1),(max(view_range)+1))
    y <- c(-0.1, order_coords[as.double(selected_node$n),3] + 1)

    x_axis$range <- x
    y_axis$range <- y

    dend_G <- graph_from_adjacency_matrix(dend_to_adj_mat(hc)) 
    es <- as.data.frame(get.edgelist(dend_G))
    
    
    dend_network <- 
      plot_ly(x = ~order_coords$x,
              y = ~order_coords$y
              ) %>% 
      add_segments(data = data.frame(
                   x = c(order_coords$x[es$V1], order_coords$x[es$V2]),
                   xend = c(order_coords$x[es$V2],order_coords$x[es$V2]), 
                   y = c(order_coords$y[es$V1], order_coords$y[es$V1]), 
                   yend = c(order_coords$y[es$V1], order_coords$y[es$V2])),
                   x = ~x,xend = ~xend,y = ~y,yend = ~yend,
                   mode='lines',color = I('black'),size = I(1), alpha = 0.5)%>%
    add_trace(type='scatter',
              mode = "markers",
              text = ~R$labels,
              marker = list(color = R$colors))
    
    dend_network<-plotly::layout(
      dend_network,
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
    net = plotly_empty() %>% add_annotations(text = "Select a Cluster to View", showarrow = F) %>% layout(xaxis = axis, yaxis = axis)
  }
  else{
      net<-plotly::layout(
      cluster_net(R, as.double(selected_node$n)),
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
        add_annotations(text = "Select a Cluster for Annotation Details", showarrow = F) %>% 
        layout(xaxis = axis, yaxis = axis)

    }

    else{
      barplot_df <- get_anno_data(R, selected_node$n)
      sunburst <- plot_ly(barplot_df, labels = ~labels, parents = ~parents, values = ~values,type = 'sunburst', branchvalues = 'remainder')
    }
    
    sunburst %>% layout(
      title="Module Node Annotations")
    sunburst
    }
  )

  output$mytable = DT::renderDataTable({
    R$annos[R$clusts[[selected_node$n]],]
  })

}



shinyApp(ui, server)


