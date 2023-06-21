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

  output$all_modules_table <- DT::renderDataTable({
    peak_df <- R$clust_info[(color_list$colors %in% c("#E7298A", "#A6D854", "#FC8D62")&R$clust_info$Size>1),] %>%
      dplyr::select(tidyselect::any_of(c("Size","pheno1_densities","pheno2_densities","densities")))
    DT::datatable(peak_df, selection="single")
  })

  tableProxy <- DT::dataTableProxy("all_modules_table")

  output$analyte_table <- DT::renderDataTable(
    data.frame(R$annos),selection="single"
  )

  analyteProxy <- DT::dataTableProxy("analyte_table")

  ### Peak Annoation Table ###
  output$anno_table <- reactable::renderReactable({
    anno_list <- colnames(R$annos)
    peak_df <- R$clust_info[(color_list$colors %in% c("#E7298A", "#A6D854", "#FC8D62")&R$clust_info$Size>1),] %>%
      dplyr::select(tidyselect::any_of(c("ClusterID","Size","pheno1_densities","pheno2_densities","densities")))
    peak_df[colnames(R$annos)] <- peak_df$ClusterID

    anno_df <<- lapply(1:nrow(peak_df), function(i){
      annos <- R$annos[R$clusts[[peak_df$ClusterID[i]]],]
      annos$cluster <- rep(peak_df$ClusterID[i], nrow(annos))
      annos %>% data.frame() %>%
        replace(is.na(.), "No Data")
    }) %>% dplyr::bind_rows()

    pal_list <<- lapply(anno_list, function(x){
      anno_vec = anno_df %>% ungroup() %>% dplyr::select(!!x) %>% unlist()
      pal = autofocus:::get_pie_color(anno_vec)
    })
    names(pal_list) <<- anno_list

    cd_fun_list <- lapply(anno_list, function(x){templates::tmpl(autofocus:::pie_plot_tmpl_fun, anno = x)})
    names_list <- lapply(anno_list, function(x){paste(x," Distribution")})
    cd_list <- lapply(1:length(anno_list), function(x){reactable::colDef(name=names_list[[x]],
                                                              cell=cd_fun_list[[x]])})

    if("densities" %in% colnames(R$clust_info)){
      columns_list <- list(reactable::colDef(maxWidth = 100),
                           reactable::colDef(maxWidth = 100),
                           reactable::colDef(maxWidth = 100))
    }else{
      columns_list <- list(reactable::colDef(maxWidth = 100),
                           reactable::colDef(maxWidth = 100),
                           reactable::colDef(maxWidth = 100),
                           reactable::colDef(maxWidth = 100))
    }

    full_columns_list <- c(columns_list, cd_list)
    names(full_columns_list) <- colnames(peak_df)

    reactable::reactable(peak_df,
                         columns = full_columns_list,
                         selection = "single",
                         onClick = "select")
  })

  output$single_module_table<-DT::renderDataTable({
    data.frame(R$annos[R$clusts[[selected_node$n]],])
  })


  ### Network/Driver section ###
  output$network <- networkD3::renderForceNetwork({
    if(!is.na(selected_node$n)){
      graph = R$graphs[[selected_node$n]]
    }else{
      graph = NA
    }
    if(length(graph)!=1){
      ColourScale_nodes <- paste0('d3.scaleOrdinal()
                          .domain(["phenotype", "analyte","confounder"])
                          .range(["',paste(unlist(color_palette[c(1,3,9)]), collapse='","'),'"]);')

      d3net<- networkD3::igraph_to_networkD3(graph, group=igraph::V(graph)$NodeType)
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
    }
  })

  output$drivers<-DT::renderDataTable(data.frame(R$annos)[rownames(R$annos)%in%autofocus:::get_drivers(R,selected_node$n),], caption="Driver information")


  # Click on a peak in the peak list
  selected_node = shiny::reactiveValues(n = NA, last=NA)
  shiny::observeEvent(reactable::getReactableState("anno_table"),{
    selected_row = reactable::getReactableState("anno_table")$selected
    if(!is.null(selected_row)){
      dt <- R$clust_info[color_list$colors %in% c("#E7298A", "#A6D854", "#FC8D62"),]
      selected_node$last = selected_node$n
      selected_node$n = dt[selected_row,1]
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
    }
  })

  # Click on an analyte in the analyte list
  selected_analyte = shiny::reactiveValues(n = NA, last=NA)
  shiny::observeEvent(input$analyte_table_row_last_clicked,{
    analyte <- input$analyte_table_row_last_clicked
    ancestors <- autofocus:::get_ancestors(R, -1*analyte, c())
    peak <- ancestors[which(color_list$colors[ancestors] %in% c("#E7298A", "#A6D854", "#FC8D62"))]
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
