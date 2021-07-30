source("functions.R")
CoCena_network_generation_UI<-function(id){
  
  ns<-NS(id)
  tagList(
    br(),
    fluidRow(column(12,DT::dataTableOutput(ns("cluster_info")))),
    br(),
    br(),
    fluidRow(column(9,plotOutput(ns("heatmap_cluster"))),
             column(3,plotOutput(ns("network_vis"))))
    
  )
}

CoCena_network_generation<-function(input,output,session,filt_cutoff,GFC_all_genes,min_nodes_cluster,voi){
  
  edgelist<-reactive({
    
    filtered_cutoff<-filt_cutoff()[["filt_cutoff_data"]]
    filtered_cutoff$weight <- filt_cutoff()[["filt_cutoff_data"]]$rval
    filtered_cutoff$rval <- NULL
    filtered_cutoff$pval <- NULL
    
    merged_net <- igraph::graph_from_data_frame(filtered_cutoff, directed=FALSE)
    merged_net <- igraph::simplify(merged_net, edge.attr.comb=list(weight="min", "ignore"))
    
    merged_net
    })
  
  cluster_calc<-reactive({
    
    withProgress(message = 'Processing:',
                 detail = 'Calculating network clusters. This may take a while...',value = 0, {
    #set prerequisits
    GFC<-GFC_all_genes()
    max_cluster_count_per_gene = 10
    min_cluster_size = as.numeric(min_nodes_cluster())
    
    
    comps = count_components(edgelist())
    
    cluster_algo="auto"
    cluster_algo_list =c("cluster_label_prop",
                         "cluster_fast_greedy",
                         "cluster_louvain",
                         "cluster_infomap",
                         "cluster_walktrap")
    
    algos_to_use = switch(cluster_algo=="auto", cluster_algo_list, "cluster_louvain")
    
    # new (actually old) seeding method:
    set.seed(168575)
    df_modularity_score = do.call("rbind", lapply(algos_to_use,
                                                  cluster_calculations,
                                                  graph_obj=edgelist(),
                                                  case="test", iter=1))
    
    
    cluster_algo_used= df_modularity_score %>%
      dplyr::filter(modularity_score==max(modularity_score)) %>%
      dplyr::select(cluster_algorithm) %>%
      as.character()
    
    igraph_list = list()
    igraph_list[[1]] = edgelist()
    
    ###apply the best clustering algorithm
    gene_which_cluster=do.call("cbind", lapply(1:20,
                                               cluster_calculations,
                                               algo=cluster_algo_used,
                                               case="best",
                                               graph_obj=edgelist()))
    
    ##frequency and identity of cluster assingment of genes
    if(base::ncol(gene_which_cluster) > 1) {
      gene_cluster_ident = apply(gene_which_cluster,1, function(x){
        if(length(unique(x)) > max_cluster_count_per_gene) {   #LISA: was >=
          0
        }else{
          names(which(table(x) == max(table(x))))[1]
        }
      })
    } else{ gene_cluster_ident = gene_which_cluster[,1]}
    
    
    
    white_genes_clustercounts <- as.integer(grep(gene_cluster_ident, pattern = "\\b0\\b") %>%
                                              length() %>% as.character())
    
    cluster_Data = data.frame(genes=vertex_attr(edgelist(), "name"),
                              clusters= paste0("Cluster ",gene_cluster_ident),
                              stringsAsFactors = FALSE)
    
    #summarize the data
    #produces a table where col are cluster name, number of components,
    #names of genes in cluster
    
    dfk=cluster_Data %>%
      dplyr::count(clusters,genes) %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(gene_no= sum(n), gene_n = paste0(genes,collapse = ",")) %>%
      dplyr::mutate(cluster_included=ifelse(gene_no>=min_cluster_size, "yes", "no"), color="white")

    ##ggplot
    color.cluster <- c("orchid", "maroon", "darkgreen",  "darkorange", "darkgrey", "gold", "steelblue", "indianred",
                       "pink", "lightgreen", "lightblue","sandybrown",   "khaki",  "turquoise","darkblue",
                       "cadetblue","greenyellow","cyan", "thistle", "darkmagenta", "coral", "red", "blue",
                       "green", "yellow", "brown", "black", "darkgoldenrod", "cornsilk", "firebrick", "deeppink",
                       "dodgerblue", "lightpink", "midnightblue", "slategray", "aquamarine",
                       "chocolate", "darkred", "navy", "olivedrab", "peachpuff",
                       "seagreen", "plum", "tomato", "snow", "wheat")
    
    plot_clust <- ggplot_build(ggplot(data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ],
                            aes(x = clusters)) +
      geom_bar(aes(fill = clusters)) +
      scale_fill_manual(values = color.cluster))
    
    #plot_clust <- ggplot_build(plot_clusters)
    
    dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color" ] <- plot_clust$data[[1]]["fill"]
    
    
    
    
    white_genes_clustersize <- as.integer(dfk %>% dplyr::filter(cluster_included=="no")%>%
                                            dplyr::summarise(n=sum(gene_no)) %>% purrr::map(1))

    # for each cluster produces a row of means per condition (info data) for all genes within the cluster

    cluster_df=dfk
    gfc_dat = GFC
    
    dfk_allinfo=do.call("rbind", lapply(1:nrow(dfk), gfc_mean_clustergene, cluster_df= dfk,
                                        gfc_dat=GFC))
    dfk_allinfo$vertexsize = ifelse(dfk_allinfo$cluster_included=="yes",3,1)
    #dfk_allinfo=as.data.frame(dfk_allinfo)
    print("dfk_allinfo")
    print(dfk_allinfo)
    dfk_allinfo
    })
  })
  
  output$cluster_info<-DT::renderDataTable({
    cluster_information<-cluster_calc()[order(cluster_calc()$cluster_included, decreasing = T),]
    DT::datatable(data=as.matrix(cluster_information),class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(columnDefs = list(list(
                    targets = 3,
                    render = JS(
                      "function(data, type, row, meta) {",
                      "return type === 'display' && data.length > 40 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
                      "}")
                  )),searchHighlight = TRUE,deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download cluster information')#I('colvis')
                                 )))
  })
  
  heatmap_cluster_reactive<-reactive({

      cluster_info=cluster_calc()
      GFCs=GFC_all_genes()
      group=voi()

      gc()
      # filter for included clusters (non-white)
      c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
      mat_heatmap <- NULL

        for (c in unique(c_df$color)){
          #get genes from the original cluster
          genes <- c_df[c_df$color == c, ] %>%
            dplyr::pull(., "gene_n") %>%
            base::strsplit(., split = ",") %>%
            unlist(.)
          
          
          # GFCs of new data set, where genes are found in original cluster
          c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
          
          if(is.vector(c_GFCs)){
            c_GFC_means <- cGFCs
          }else{
            c_GFC_means <- apply(c_GFCs[, c(1:(ncol(c_GFCs)-1))]%>% as.data.frame(), 2, base::mean)
          }
          
          mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
          
        }
        rownames(mat_heatmap) <- c_df$color
    
      colnames(mat_heatmap) <- colnames(GFCs)[1:(ncol(GFCs)-1)]

        cluster_colors <- factor(c_df$color)
        names(cluster_colors) <- c_df$color
        row_order <- unique(c_df$color)

        ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors,
                                                                      simple_anno_size = unit(0.5, "cm")),
                                                genes = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                                gene_nums = anno_text(c_df$gene_no, width = unit(1.5, "cm"), gp = gpar(fontsize = 12)),
                                                
                                                
                                                which = "row",
                                                width = unit(4.5, "cm"),
                                                annotation_name_side = "top",
                                                gap = unit(2, "mm"),
                                                annotation_name_rot = 0,
                                                annotation_name_gp = gpar(fontsize = 12))
        
        lgd_list <- list()
      
      anno_list <- NULL

      all_conditions <- NULL
      
      hm <- ComplexHeatmap::Heatmap(mat_heatmap,
                                    right_annotation = ha,
                                    col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1))),
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_columns = "euclidean",
                                    clustering_method_rows = "complete",
                                    clustering_method_columns = "complete",
                                    cluster_columns = T,
                                    cluster_rows = T,
                                    column_names_rot = 90,
                                    column_names_centered = F,
                                    row_names_gp = gpar(fontsize = 12),
                                    column_names_gp = gpar(fontsize = 12),
                                    rect_gp = grid::gpar(col = "black"),
                                    heatmap_legend_param = list(title = "", legend_height = unit(3, "cm")), column_km = 0)
      
      hm_w_lgd <- ComplexHeatmap::draw(hm %v% anno_list, annotation_legend_list = lgd_list, merge_legends = T,
                                       padding = unit(c(2, 2, 2, 30), "mm"))
      
      hm_w_lgd
      
  })
  
  output$heatmap_cluster<-renderPlot({
    heatmap_cluster_reactive()
  })
  
  calc_network<-reactive({
    output<-list()
    
      gene_to_cluster <- do.call(rbind, apply(cluster_calc(),1, function(x){
        tmp <- x["gene_n"] %>%
          base::strsplit(., split = ",")%>%
          unlist(.)
        data.frame(name = tmp, color = rep(x["color"], length(tmp)))
      }))
      gene_to_cluster$name <- as.character(gene_to_cluster$name)
      gene_to_cluster$color <- as.character(gene_to_cluster$color)
      
      network=edgelist()
      
      if(length(V(network)$name[!V(network)$name %in% gene_to_cluster$name]) > 0){
        gene_to_cluster <- rbind(gene_to_cluster, data.frame(name = V(network)$name[!V(network)$name %in% gene_to_cluster$name],
                                                             color = "white"))
        gene_to_cluster$name <- as.character(gene_to_cluster$name)
        gene_to_cluster$color <- as.character(gene_to_cluster$color)
      }
      gene_to_cluster <- gene_to_cluster[match(V(network)$name, gene_to_cluster$name),]
      V(network)$color <- gene_to_cluster$color
      if("white" %in% unique(gene_to_cluster$color)){
        network <- delete.vertices(network, gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))
      }


        set.seed(123)
        l <- igraph::layout.fruchterman.reingold(network) 
        rownames(l) <- V(network)$name
        #l <- "layout_with_fr"
        
        output[["labelled_network"]] <- list()
        output[["network_col_by_module"]] <- network

        output[["plotted_network"]]<-igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                            layout = l, main = "co-expression network coloured by cluster")
        
        output
      
  })
  
  output$network_vis<-renderPlot({
    calc_network()[["plotted_network"]]
  })
  
  return(list(cluster_calc=reactive({cluster_calc()})))
}