suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dendextend))
source(codes.makepath("autofocus/Frontend_Modules.R"))


# ui <- fluidPage(theme = shinytheme("journal"),
#   titlePael("AutoFocus Run Results"),
#       fluidRow(
#         column(12,plotlyOutput("dendro", height="50%")),
#                #plotlyOutput("barcharts")),
#         column(12,dataTableOutput("all_modules_table",height="50%"),
#                dataTableOutput("single_module_table",height="50%"))
# ))
load("~/Box/results/Annalise/QD_R.rda")
R<-QD_results_lm

body <- dashboardBody(
  fluidRow(tabBox(title=NULL,width=12,
                  tabPanel("Tree View",plotlyOutput("dendro")),
                  tabPanel("List View",dataTableOutput("all_modules_table")))),
    fluidRow(
            dataTableOutput("single_module_table",height="50%")
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
                                          marker = list(color=c(R$colors[selected_node$last],"yellow"), size = 9)))
    
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
                                                       marker = list(color=c(R$colors[selected_node$last],"yellow"), size = 9)))
    }
  })
    
  
  ### Annotation section ###
  # output$barplots <- renderPlotly({
  #   axis <- list(
  #     
  #     showline = TRUE,
  #     showticklabels = F,
  #     mirror="ticks",
  #     showgrid = F,
  #     zeroline = F,
  #     linecolor = toRGB("black"),
  #     linewidth = 2
  #   )
  #   if (is.na(selected_node$n)){
  #     sunburst = plotly_empty() %>% 
  #       add_annotations(text = "Select a Cluster for Annotation Details", showarrow = F) %>% 
  #       layout(xaxis = axis, yaxis = axis)
  # 
  #   }
  # 
  #   else{
  #     barplot_df <- get_anno_data(R, selected_node$n)
  #     sunburst <- plot_ly(barplot_df, labels = ~labels, parents = ~parents, values = ~values,type = 'sunburst', branchvalues = 'remainder')
  #   }
  #   
  #   sunburst %>% layout(
  #     title="Module Node Annotations")
  #   sunburst
  #   }
  # )


}

shinyApp(ui, server)

