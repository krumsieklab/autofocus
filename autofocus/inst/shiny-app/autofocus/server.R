suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(magrittr))

server <- function(input, output) {

  y_axis <- list(
    title = "",
    showline = TRUE,
    showticklabels = F,
    mirror="ticks",
    showgrid = F,
    zeroline = F,
    linecolor = plotly::toRGB("black"),
    linewidth = 2
  )
  x_axis <- list(
    title = "",
    showline = TRUE,
    showticklabels = F,
    mirror="ticks",
    showgrid = F,
    zeroline = F,
    linecolor = plotly::toRGB("black"),
    linewidth = 2
  )

  color_list<-shiny::reactiveValues(colors = NA)
  shiny::observe({
    color_list$colors =  autofocus:::peak_finder_wrapper(R, input$threshold)
  })

  output$dendro <- plotly::renderPlotly({

    ### Dendrogram ###
    dend_network<-plotly::layout(
      autofocus:::plot_dend(R, input$threshold),
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = F
    )
    dend_network$source = "A"
    dend_network
  })

  dendroProxy <- plotly::plotlyProxy("dendro")

  output$all_modules_table <- DT::renderDataTable(
    R$clust_info[color_list$colors=="red",], selection="single"
  )

  tableProxy <- DT::dataTableProxy("all_modules_table")

  output$analyte_table <- DT::renderDataTable(
    data.frame(R$annos),selection="single"
  )

  analyteProxy <- DT::dataTableProxy("analyte_table")


  output$single_module_table<-DT::renderDataTable(data.frame(R$annos[R$clusts[[selected_node$n]],]))


  ### Annotation section ###
  output$barplots <- networkD3::renderSankeyNetwork({
    mat <- R$annos[R$clusts[[selected_node$n]],]
    print(dim(mat))
    if(dim(mat)[2]==1){
      mat <- cbind(mat, mat)
    }
    df_mat<- mat %>% data.frame
    for(n in 1:ncol(df_mat)) {                                   # Replace NA in all columns
      df_mat[ , n][is.na(df_mat[ , n])] <- paste0("Missing ",colnames(df_mat)[n])
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
    sn <- networkD3::sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                        Target = 'target', Value = 'count', NodeID = 'name', fontSize=10)

    # add origin back into the links data because sankeyNetwork strips it out
    sn$x$links$origin <- links$origin


    # add onRender JavaScript to set the click behavior
    sn<-htmlwidgets::onRender(
      sn,
      paste0('
      function(el, x) {
        var nodes = d3.selectAll(".node");
        var links = d3.selectAll(".link");
        var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i).sort(function(a, b){return a - b});
        var labels = ["',paste(unlist(colnames(df_mat)), collapse='","'),'"];
        cols_x.forEach((d, i) => {
          d3.select(el).select("svg")
            .append("text")
            .attr("x", d)
            .attr("y", 12)
            .text(labels[i]);
        })
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
      ))
    sn
  })

  ### Network/Driver section ###
  output$network <- networkD3::renderForceNetwork({
    ColourScale_nodes <- paste0('d3.scaleOrdinal()
                          .domain(["phenotype", "analyte","confounder"])
                          .range(["',paste(unlist(color_palette[c(1,3,9)]), collapse='","'),'"]);')

    graph = R$graphs[[selected_node$n]]
    d3net<- networkD3::igraph_to_networkD3(graph, group=V(graph)$NodeType)
    drivers <- igraph::V(graph)$Driver
    networkD3::forceNetwork(d3net$links,
                 d3net$nodes,
                 NodeID="name",
                 Group="group",
                 colourScale=networkD3::JS(ColourScale_nodes),
                 linkColour = ifelse(rowSums(d3net$links==0), color_palette[4],color_palette[10]),
                 zoom=T,
                 opacity=1,
                 legend=T,
    )
  })

  output$drivers<-DT::renderDataTable(data.frame(R$annos)[rownames(R$annos)%in%autofocus:::get_drivers(R,selected_node$n),], caption="Driver information")


  # Click on a peak in the peak list
  selected_node = shiny::reactiveValues(n = NA, last=NA)
  shiny::observeEvent(input$all_modules_table_row_last_clicked,{
    dt <- R$clust_info[color_list$colors=="red",]
    selected_node$last = selected_node$n
    selected_node$n = dt[input$all_modules_table_row_last_clicked,1]
    old_label <- if(is.na(selected_node$last)) "" else autofocus:::get_node_label(R,selected_node$last)
    new_label <- autofocus:::get_node_label(R, selected_node$n)
    view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
    x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
    y <- R$clust_info[selected_node$n,]$Coord_Y
    y_axis$range <- c(-(y*0.1), y+(y*0.1))
    plotly::plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>%
      plotly::plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, R$clust_info[selected_node$n,]$Coord_X),
                                          y=c(R$clust_info[selected_node$last,]$Coord_Y, R$clust_info[selected_node$n,]$Coord_Y),
                                          type="scatter",
                                          mode="markers",
                                          text = c(old_label,new_label),
                                          marker = list(color=c(color_list$colors[selected_node$last],"yellow"), size = 9)))

  })

  # Click on an analyte in the analyte list
  selected_analyte = shiny::reactiveValues(n = NA, last=NA)
  shiny::observeEvent(input$analyte_table_row_last_clicked,{
    analyte <- input$analyte_table_row_last_clicked
    ancestors <- autofocus:::get_ancestors(R, -1*analyte, c())
    peak <- ancestors[which(color_list$colors[ancestors]=="red")]
    if(length(peak)==0){
      shinyalert::shinyalert(
        title = "Insignificant", type = "warning",
        text="This analyte does not belong to a significant module, please choose another"
      )
    }
    else{
      selected_analyte$last = selected_analyte$n
      selected_analyte$n = input$analyte_table_row_last_clicked+dim(R$HCL$merge)[1]
      old_analyte_label <- if(is.na(selected_analyte$last)) "" else autofocus:::get_node_label(R,selected_analyte$last)
      new_analyte_label <- autofocus:::get_node_label(R, selected_analyte$n)
      selected_node$last = selected_node$n
      selected_node$n=peak
      old_label <- if(is.na(selected_node$last)) "" else autofocus:::get_node_label(R,selected_node$last)
      new_label <- autofocus:::get_node_label(R, selected_node$n)
      view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
      x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
      y <- R$clust_info[selected_node$n,]$Coord_Y
      y_axis$range <- c(-(y*0.1), y+(y*0.1))
      plotly::plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>%
        plotly::plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X,
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
  shiny::observeEvent(plotly::event_data(event = "plotly_click", source="A"),{
    d<-plotly::event_data("plotly_click", source="A")
    if (d$y != 0){
      selected_node$last = selected_node$n
      selected_node$n = which(R$clust_info$Coord_Y ==round(d$y,digits=10))
      old_label <- if(is.na(selected_node$last)) "" else autofocus:::get_node_label(R,selected_node$last)
      new_label <- autofocus:::get_node_label(R, selected_node$n)
      print(c(old_label,new_label))
      DT::selectRows(tableProxy, selected_node$n)
      view_range <- match(R$clusts[as.double(selected_node$n)][[1]], R$HCL$order)
      x_axis$range <- c((min(view_range)-1),(max(view_range)+1))
      y <- R$clust_info[selected_node$n,]$Coord_Y
      y_axis$range <- c(-(y*0.1), y+(y*0.1))
      plotly::plotlyProxyInvoke(dendroProxy, "relayout", list(xaxis=x_axis, yaxis=y_axis)) %>%
        plotly::plotlyProxyInvoke("addTraces", list(x=c(R$clust_info[selected_node$last,]$Coord_X, R$clust_info[selected_node$n,]$Coord_X),
                                            y=c(R$clust_info[selected_node$last,]$Coord_Y, R$clust_info[selected_node$n,]$Coord_Y),
                                            type="scatter",
                                            mode="markers",
                                            text = c(old_label, new_label),
                                            marker = list(color=c(color_list$colors[selected_node$last],"yellow"))))

    }
  })

}
