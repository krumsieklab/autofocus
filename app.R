source(codes.makepath("autofocus/Frontend_Modules.R"))


args = commandArgs(trailingOnly = T)
load(file=args[1])

columns_to_include <- args[2:length(args)]
R$annos<-R$annos[,columns_to_include]
body <- dashboardBody(
  fluidRow(tabBox(title=NULL,width=12,
                  tabPanel("Tree View",plotlyOutput("dendro")),
                  tabPanel("List View",dataTableOutput("all_modules_table")))),
    fluidRow(
            dataTableOutput("single_module_table"),
            plotlyOutput("barplots")
  )
)

ui<-dashboardPage(dashboardHeader(title="AutoFocus Run Results"), 
                  dashboardSidebar(disable = TRUE), 
                  body)


server <- function(input, output) {

  output$dendro <- renderPlotly({

    ### Dendrogram ###
    dend_network<-plotly::layout(
      plot_dend(R),
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = F
    )
    dend_network
   
  })
  
  output$all_modules_table <- DT::renderDataTable(
    data.frame(R$clust_info)[R$clust_info$Size>1,], selection="single"
  )
  
  #output$barcharts<-renderPlotly()
  
  output$single_module_table<-DT::renderDataTable(data.frame(R$annos[R$clusts[[selected_node$n]],]))
  
  dendroProxy <- plotlyProxy("dendro")
  
  tableProxy = dataTableProxy("all_modules_table")
  
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
  
  selected_node = reactiveValues(n = NA, last=NA)
  observeEvent(input$all_modules_table_row_last_clicked,{
    selected_node$last = selected_node$n
    selected_node$n = input$all_modules_table_row_last_clicked
    view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
    x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
    y <- R$clust_info[selected_node$n,]$Coord_Y
    y_axis$range <- c(-(y*0.1), y+(y*0.1))
    plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>% 
      plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, R$clust_info[selected_node$n,]$Coord_X),
                                          y=c(R$clust_info[selected_node$last,]$Coord_Y, R$clust_info[selected_node$n,]$Coord_Y),
                                          type="scatter",
                                          mode="markers",
                                          marker = list(color=c(R$clust_info$colors[selected_node$last],"yellow"), size = 9)))
    
  })
  observeEvent(event_data("plotly_click"),{
    d<-event_data("plotly_click")
    if (d$y != 0){
      selected_node$last = selected_node$n
      selected_node$n = which(R$clust_info$Coord_Y ==round(d$y,digits=10))
      selectRows(tableProxy, selected_node$n)
      view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
      x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
      y <- R$clust_info[selected_node$n,]$Coord_Y
      y_axis$range <- c(-(y*0.1), y+(y*0.1))
      plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>% 
      plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, R$clust_info[selected_node$n,]$Coord_X),
                                                       y=c(R$clust_info[selected_node$last,]$Coord_Y, R$clust_info[selected_node$n,]$Coord_Y),
                                                       type="scatter",
                                                       mode="markers",
                                                       marker = list(color=c(R$clust_info$colors[selected_node$last],"yellow"), size = 9)))
    }
  })
    
  
  ### Annotation section ###
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
      barplots <- get_anno_data(R, selected_node$n)
      sunburst <-subplot(barplots, margin=0.09)
    }

    sunburst %>% layout(
      title="Module Node Annotations")
    sunburst
    }
  )


}

shinyApp(ui, server)

