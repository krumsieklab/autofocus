source(codes.makepath("autofocus/Frontend_Modules.R"))
library(networkD3)
library(htmlwidgets)
suppressPackageStartupMessages(library(MetBrewer))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinyalert))

#args = commandArgs(trailingOnly = T)
#load(file=args[1])

#columns_to_include <- args[2:length(args)]
#R$annos<-R$annos[,columns_to_include]

R$annos <- R$annos[,c("platform","SUPER_PATHWAY","SUB_PATHWAY")]
color_palette<-met.brewer("Hiroshige",10)
body <- dashboardBody(
  fluidRow(tabBox(title=NULL,
                  width=12,
                  tabPanel("Tree View",plotlyOutput("dendro")),
                  tabPanel("Peak List",dataTableOutput("all_modules_table")),
                  tabPanel("Analyte List",dataTableOutput("analyte_table")),
                  sliderInput("threshold", "Density Threshold:", min=0, max=1, value=0.95)),
          tabBox(title = NULL, 
                  width=12,
                  tabPanel("Module members", dataTableOutput("single_module_table")),
                  tabPanel("Module Network and Drivers", 
                    fluidRow(column(forceNetworkOutput("network"),width=9), column(dataTableOutput("drivers"),width=11))),
                  tabPanel("Annotation plots", sankeyNetworkOutput("barplots"))))
)

ui<-dashboardPage(dashboardHeader(title="AutoFocus Run Results"), 
                  dashboardSidebar(disable = TRUE), 
                  body)


server <- function(input, output) {

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
  
  color_list<-reactiveValues(colors = NA)
  observe({
    color_list$colors =  peak_finder_wrapper(R, input$threshold)
  })
  
  output$dendro <- renderPlotly({

    ### Dendrogram ###
    dend_network<-plotly::layout(
      plot_dend(R, input$threshold),
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = F
    )
    dend_network$source = "A"
    dend_network
  })
  
  dendroProxy <- plotlyProxy("dendro")
  
  output$all_modules_table <- DT::renderDataTable(
    R$clust_info[color_list$colors=="red",], selection="single"
  )
  
  tableProxy <- dataTableProxy("all_modules_table")
  
  output$analyte_table <- DT::renderDataTable(
    data.frame(R$annos),selection="single"
  )
  
  analyteProxy <- dataTableProxy("analyte_table")
  
  
  output$single_module_table<-DT::renderDataTable(data.frame(R$annos[R$clusts[[selected_node$n]],]))
  
  
  ### Annotation section ###
  output$barplots <- renderSankeyNetwork({
      mat <- R$annos[R$clusts[[selected_node$n]],]
      df_mat<- mat %>% data.frame
      for(i in 1:ncol(df_mat)) {                                   # Replace NA in all columns
        df_mat[ , i][is.na(df_mat[ , i])] <- paste0("Missing ",colnames(df_mat)[i])
      }
      links <-
        df_mat %>%
        unite(origin, colnames(df_mat),remove=F) %>%
        mutate(row = row_number()) %>%
        gather("column", "source", -row, -origin) %>%
        mutate(column = match(column, names(df))) %>%
        arrange(row, column) %>%
        group_by(row) %>%
        mutate(target = lead(source)) %>%
        ungroup() %>%
        filter(!is.na(target)) %>%
        select(source, target, origin) %>%
        group_by(source, target, origin) %>%
        summarise(count = n()) %>%
        ungroup()
      
      nodes <- data.frame(name = unique(c(links$source, links$target)))
      links$source <- match(links$source, nodes$name) - 1
      links$target <- match(links$target, nodes$name) - 1
      sn <- sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                          Target = 'target', Value = 'count', NodeID = 'name', fontSize=10)
      
      # add origin back into the links data because sankeyNetwork strips it out
      sn$x$links$origin <- links$origin
      
      
      # add onRender JavaScript to set the click behavior
      sn<-htmlwidgets::onRender(
        sn,
        '
      function(el, x) {
        var nodes = d3.selectAll(".node");
        var links = d3.selectAll(".link");
        nodes.on("mousedown.drag", null); // remove the drag because it conflicts
        nodes.on("click", clicked);
        function clicked(d, i) {
          links
            .style("stroke-opacity", function(d1) {
                var str = d1.origin
                return str.includes(d.name) ? 0.5 : 0.2;
              });
        }
      }
      '
      )
      sn
  })
  
  ### Network/Driver section ###
  output$network <- renderForceNetwork({
    ColourScale_nodes <- paste0('d3.scaleOrdinal()
                          .domain(["phenotype", "analyte","confounder"])
                          .range(["',paste(unlist(color_palette[c(1,3,9)]), collapse='","'),'"]);')
    
    graph = R$graphs[[selected_node$n]]
    d3net<- igraph_to_networkD3(graph, group=V(graph)$NodeType)
    drivers <- V(graph)$Driver
    forceNetwork(d3net$links, 
                 d3net$nodes, 
                 NodeID="name",
                 Group="group", 
                 colourScale=JS(ColourScale_nodes),
                 linkColour = ifelse(rowSums(d3net$links==0), color_palette[4],color_palette[10]),
                 zoom=T,
                 opacity=1,
                 legend=T,
                 )
  })

  output$drivers<-DT::renderDataTable(data.frame(R$annos)[rownames(R$annos)%in%get_drivers(R,selected_node$n),], caption="Driver information")
  

  # Click on a peak in the peak list
  selected_node = reactiveValues(n = NA, last=NA)
  observeEvent(input$all_modules_table_row_last_clicked,{
    dt <- R$clust_info[color_list$colors=="red",]
    selected_node$last = selected_node$n
    selected_node$n = dt[input$all_modules_table_row_last_clicked,1]
    old_label <- if(is.na(selected_node$last)) "" else get_node_label(R,selected_node$last)
    new_label <- get_node_label(R, selected_node$n)
    view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
    x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
    y <- R$clust_info[selected_node$n,]$Coord_Y
    y_axis$range <- c(-(y*0.1), y+(y*0.1))
    plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>%
      plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, R$clust_info[selected_node$n,]$Coord_X),
                                          y=c(R$clust_info[selected_node$last,]$Coord_Y, R$clust_info[selected_node$n,]$Coord_Y),
                                          type="scatter",
                                          mode="markers",
                                          text = c(old_label,new_label),
                                          marker = list(color=c(color_list$colors[selected_node$last],"yellow"), size = 9)))

  })
  
  # Click on an analyte in the analyte list
  selected_analyte = reactiveValues(n = NA, last=NA)
  observeEvent(input$analyte_table_row_last_clicked,{
    analyte <- input$analyte_table_row_last_clicked
    ancestors <- get_ancestors(R, -1*analyte, c())
    peak <- ancestors[which(color_list$colors[ancestors]=="red")]
    if(length(peak)==0){
      shinyalert(
        title = "Insignificant", type = "warning",
        text="This analyte does not belong to a significant module, please choose another"
        )
    }
    else{
      selected_analyte$last = selected_analyte$n
      selected_analyte$n = input$analyte_table_row_last_clicked+dim(R$HCL$merge)[1]
      old_analyte_label <- if(is.na(selected_analyte$last)) "" else get_node_label(R,selected_analyte$last)
      new_analyte_label <- get_node_label(R, selected_analyte$n)
      selected_node$last = selected_node$n
      selected_node$n=peak
      old_label <- if(is.na(selected_node$last)) "" else get_node_label(R,selected_node$last)
      new_label <- get_node_label(R, selected_node$n)
      view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
      x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
      y <- R$clust_info[selected_node$n,]$Coord_Y
      y_axis$range <- c(-(y*0.1), y+(y*0.1))
      plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>%
        plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, 
                                                R$clust_info[selected_node$n,]$Coord_X,
                                                R$clust_info[selected_analyte$last,]$Coord_X,
                                                R$clust_info[selected_analyte$n,]$Coord_X),
                                            y=c(R$clust_info[selected_node$last,]$Coord_Y, 
                                                R$clust_info[selected_node$n,]$Coord_Y,
                                                R$clust_info[selected_analyte$last,]$Coord_Y,
                                                R$clust_info[selected_analyte$n,]$Coord_Y),
                                            type="scatter",
                                            mode="markers",
                                            text = c(old_label,new_label, old_analyte_label, new_analyte_label),
                                            marker = list(color=c(color_list$colors[selected_node$last],
                                                                  "yellow",
                                                                  color_list$colors[selected_analyte$last],
                                                                  "blue"))))
      
      }
    })
  
  # Click on a node in the dendrogram
  observeEvent(event_data(event = "plotly_click", source="A"),{
    d<-event_data("plotly_click", source="A")
    if (d$y != 0){
      selected_node$last = selected_node$n
      selected_node$n = which(R$clust_info$Coord_Y ==round(d$y,digits=10))
      old_label <- if(is.na(selected_node$last)) "" else get_node_label(R,selected_node$last)
      new_label <- get_node_label(R, selected_node$n)
      print(c(old_label,new_label))
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
                                            text = c(old_label, new_label),
                                            marker = list(color=c(color_list$colors[selected_node$last],"yellow"))))
      
    }
  })
  
}

shinyApp(ui, server)

