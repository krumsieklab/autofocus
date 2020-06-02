suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
source(codes.makepath("autofocus/internal_funcs.R"))
ui <- fluidPage(
  titlePanel("AutoFocus Results"),

  selectInput(inputId = "outcome", label = "Select Disease Phenotype",
              choices = c("No outcome", phenotypes),
              selected = "No outcome"), #Outcome

      fluidRow(mainPanel("Module Browser", plotlyOutput("dendro"), plotlyOutput("clust"))),
      fluidRow( plotlyOutput("barplots")),DT::dataTableOutput("mytable")
)


server <- function(input, output) {
  if (length(phenotypes)>1){
    R <- R_list[[which(phenotypes == input$outcome)]]
  }
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
   
    
    dend_network<-plotly::layout(
      plot_dend(R, order_coords),
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
  
  ### sunburst plot ###
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

