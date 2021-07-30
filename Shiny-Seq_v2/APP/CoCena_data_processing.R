###data processing
source("functions.R")
CoCena_data_processing_UI <- function(id){
  
  ns<-NS(id)
  tagList(
    br(),
    fluidRow(column(12,DT::dataTableOutput(ns("data_processing1")))),
    #fluidRow(column(10,bsAlert(ns("cutoff")))),
    br(),
    fluidRow(column(1,actionButton(ns("set_cutoffs"), "Set custom cut-off"))),
    br(),
    fluidRow(column(6, plotOutput(ns("deg_dist_opt"))),
             column(6, plotOutput(ns("deg_dist_cust")))),
    br(),
    fluidRow(column(4,radioButtons(ns("select_cutoff"), "Which cut-off should be used", choices=list("optimal cut-off", "custom cut-off")))),
    bsModal(ns("final_cutoffs"),"Set the final cut-offs", ns("set_cutoffs"),size="large",
            helpText("Choose a minimal correlation cutoff based on which 
            the gene-co-expression network is constructed. Gene pairs are only considered as co-expressed if their 
            correlation exceeds the chosen cutoff. Their co-expression is then represented by an edge in the subsequent 
            network. To aid the cutoff selection, the most important stats for each cutoff are represented as the plot below. 
            Displayed are the R.squared value representing how strongly the scale-free topology is maintianed in the 
            resulting network, the number of edges and the number of nodes in the resulting network and the number of 
            graph components."),
            br(),
            fluidRow(column(12,plotlyOutput(ns("cutoff_plot")))),
            br(),
            fluidRow(column(4,textInput(ns("cutoff_custom"), "Choose the final minimal correlation", placeholder = "i.e. 0.972"))),
            br(),
            helpText("Based on the chosen cutoff, a linear regression is performed on the double logged degree 
                     distribution. The better the linear fit, the more the resulting network follows the scale-free 
                     topology criterion"),
            fluidRow(column(12, plotOutput(ns("deg_dist"))))) #plot the degree distribution
    
  )
}

CoCena_data_processing <- function(input,output,session,top_var,range_corr,min_corr,min_nodes_network,batch_choice,batch_corrected,normal){
  
  top_var_output <- reactive({
    withProgress(message = 'Processing:',
                 detail = 'Filtering data for most variable genes',value = 0, {
                  
                   output<-list()
                   if(as.numeric(batch_choice())==1) count_table<-normal()
                   else count_table<-batch_corrected()
                   
                   ds = count_table[order(apply(count_table,1, stats::var), decreasing=T),]
                   output[["ds"]]<-ds
                   
                   dd2 <- head(ds, as.numeric(top_var()))
                   output[["top_var_output"]]<-dd2
                   output
                 })
  })
  
  correlation <- reactive({
    withProgress(message = 'Processing:',
                 detail = 'Calculating gene correlations. This may take a while...',value = 0, {
                   dd2 = t(top_var_output()[["top_var_output"]])
                   gc()
                   corr_calc_out <- correlation_actions(dd2, as.numeric(range_corr()), as.numeric(min_corr()))
                   corr_calc_out
                 })
  })
  
  output$data_processing1 <- DT::renderDataTable(server=FALSE,{
    
    
    DT::datatable(top_var_output()[["top_var_output"]],class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download table of top variable genes data')#I('colvis')
                                 )))
    
  })
    
  cutoff_statistics <- reactive({ 
    withProgress(message = 'Processing:',
                 detail = 'Calculating correlation statistics. This may take a while...',value = 0, {
                   
                   output <- list()
                   cutoff_stats = do.call("rbind", lapply(X = correlation()[["range_cutoff"]],
                                         FUN = cutoff_prep,
                                         corrdf_r = correlation()[["correlation_df_filt"]],
                                         min_nodes_network = as.numeric(min_nodes_network())))
                   
                   output[["cutoff_stats"]] <- cutoff_stats
                   
                   cutoff_calc_out <- optcut_fun(cutoff_stats = cutoff_stats)
                   
                   output[["cutoff_calc_out"]] <- cutoff_calc_out
                  
                 })
            output
  })
  
  output$cutoff_plot<-renderPlotly({
    p <- plot_cutoffs_internal(cutoff_stats = cutoff_statistics()[["cutoff_calc_out"]][["cutoff_stats_concise"]],
                               hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL))
    p
  })
  
  deg_dist_custom<-reactive({
    req(input$cutoff_custom)
    cut_off_set<-as.numeric(input$cutoff_custom)
    if(!as.numeric(input$cutoff_custom) %in% cutoff_statistics()[["cutoff_stats"]]$cutoff){
      tmp_vec <- cutoff_statistics()[["cutoff_stats"]]$cutoff
      cut_off_set <- tmp_vec[which.min(abs(tmp_vec - as.numeric(input$cutoff_custom)))]
      
    }
    
    stats_calculated_optimal_cutoff <- cutoff_statistics()[["cutoff_stats"]][cutoff_statistics()[["cutoff_stats"]]$cutoff == cut_off_set, c("degree", "Probs")]
    
    stats <- cutoff_statistics()[["cutoff_calc_out"]][["cutoff_stats_concise"]]
    stats$cut_off <- rownames(stats)
    stats <- stats[stats$cut_off == cut_off_set,]
    
    dd_plot_calculated_optimal = ggplot(stats_calculated_optimal_cutoff,aes(x=log(degree), y= log(Probs))) +
      geom_point() +
      geom_smooth(method="lm") +
      theme_bw() + 
      ggtitle(paste0("Custom cut-off: ",cut_off_set, "; R.squared: ", round(stats[1],3), "; no. edges: ",
                     stats[2], "; no. nodes: ", stats[3], "; no. networks: ", stats[4]))+
      theme(plot.title = element_text(size = 10))
    
    dd_plot_calculated_optimal
  })
  
  output$deg_dist<-renderPlot({
    deg_dist_custom()
  })
  
  output$deg_dist_cust<-renderPlot({
    deg_dist_custom()
  })
  
  output$deg_dist_opt<-renderPlot({
    p<-cutoff_statistics()[["cutoff_calc_out"]][["dd_plot_calculated_optimal"]]
    p
    
  })
  
  cutoff_choice<-reactive({
    if(input$select_cutoff=="custom cut-off") as.numeric(input$cutoff_custom)
    else if(input$select_cutoff=="optimal cut-off") as.numeric(cutoff_statistics()[["cutoff_calc_out"]][["optimal_cutoff"]])
  })
  
  
return(list(top_var_output=reactive({top_var_output()}),
            gene_corr=reactive({correlation()}),
            corr_stats=reactive({cutoff_statistics()}),
            cutoff_choice=reactive({cutoff_choice()})))
  
}

