###CoCena data processing part 2
source("functions.R")
CoCena_data_processing2_UI<-function(id){
  
  ns<-NS(id)
  tagList(
    br(),
    fluidRow(column(6,plotOutput(ns("condition_heatmap"), height="800px", width="auto")),
             column(6,DT::dataTableOutput(ns("GFC_table"), height="800px", width="auto")))
  )
  
}

CoCena_data_processing2<-function(input,output,session,anno,top_var_output,cutoff_choice,gene_corr,voi,range_GFC){

  heatmap_CoCena <- reactive({  
    output <- list()
    info_dataset=anno()
    filt_cutoff_data = gene_corr()[["correlation_df_filt"]] %>% dplyr::filter(., rval > cutoff_choice())
    filt_cutoff_graph = igraph::graph_from_data_frame(filt_cutoff_data,directed=FALSE)
    filt_cutoff_counts = top_var_output()[["ds"]][row.names(top_var_output()[["ds"]]) %in% names(V(filt_cutoff_graph)),]
    corresp_info = info_dataset[rownames(t(top_var_output()[["top_var_output"]]))%in%rownames(info_dataset),]
  
    output[["filt_cutoff_graph"]] <- filt_cutoff_graph
    output[["filt_cutoff_data"]] <- filt_cutoff_data

    annotation_c <-info_dataset[voi()]
    rownames(annotation_c)<-colnames(filt_cutoff_counts)
      heatmap_filtered_counts <- pheatmap::pheatmap(mat = filt_cutoff_counts ,
                                                  color=rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                                  scale="row",
                                                  cluster_rows=T,
                                                  cluster_cols=T,
                                                 # annotation_colors = mat_colors,
                                                  annotation_col=annotation_c,
                                                  fontsize = 8,
                                                  show_rownames = F, 
                                                  show_colnames = T, annotation_names_col = T, 
                                                  clustering_distance_cols="euclidean", clustering_method="complete")

  
    output[["heatmap"]] <- heatmap_filtered_counts
    
    output
  })
  
  output$condition_heatmap<-renderPlot({
    heatmap_CoCena()[["heatmap"]]
  })
  
  GFC_genes<-reactive({
    withProgress(message = 'Processing:',
                 detail = 'Calculating group fold changes.',value = 0, {
    info_dataset<-anno()
    rownames(info_dataset)<-info_dataset[,1]
    if(intersect(voi(),colnames(info_dataset))%>% length()>0){

      
      info_dataset$grpvar =purrr::pmap(info_dataset[intersect(voi(),colnames(info_dataset))],
                                       paste, sep="-") %>% unlist()
    } else {

      info_dataset$grpvar = info_dataset[,1]
      
    }

      # GFC calculation with foldchange from mean
      norm_data_anno = merge(t(top_var_output()[["top_var_output"]]), info_dataset["grpvar"], by="row.names", all.x=T)
      norm_data_anno = norm_data_anno[,-1]
      
      norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
      
      trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])

      trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , base::mean)))
      trans_norm <- cbind(trans_norm , rowMeans(trans_norm))
      
      colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
      grplist=colnames(trans_norm)[-(ncol(trans_norm))]

      GFC_all_genes=do.call("cbind", lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[,"group_mean"],
                                            range_GFC = as.numeric(range_GFC())))

      GFC_all_genes= round(GFC_all_genes,3)
      GFC_all_genes$Gene = rownames(GFC_all_genes)
      })
    
      GFC_all_genes

  })
  
  output$GFC_table<-renderDataTable({
    
    DT::datatable(GFC_genes(),class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download table of group fold change data')#I('colvis')
                                 )))
  })
  
  return(list(filt_cutoff=reactive({heatmap_CoCena()}),
              GFC_all_genes=reactive({GFC_genes()})))
}