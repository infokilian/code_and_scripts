###global settings
CoCena_module_UI<-function(id){
  ns<-NS(id)
  tagList(
    br(),
    fluidRow(column(3,textInput(ns("min_nodes_network"), "Define the minimum number of nodes required for the network", value = 40)),
             column(3,br(),textInput(ns("top_var"), "Choose the number of top variable genes", value=5000))),
    br(),
    fluidRow(column(3,textInput(ns("min_nodes_cluster"), "Define the minimum number of nodes required to form a cluster", value = 20)),
             column(3,textInput(ns("min_correlation"), "Set the minimal correlation required for genes to be selected", value = 0.7))),
    br(),
    fluidRow(column(3,br(),uiOutput(ns("variable_of_interest"))),
             column(3,textInput(ns("range_correlation"), "How many correlation cut-offs should be tested", value=400))),
    br(),
    fluidRow(column(3,sliderInput(ns("range_gfc"), "Define the range of the group fold change (GFC)", min=1, max=10, value=2, step=0.1))),

  )
}

CoCena_module<-function(input,output,session,anno){
  
  output$variable_of_interest <- renderUI({
    variable_anno<-list()
    for(i in 1:(ncol(anno())-1)){
      variable_anno[[i]]<-colnames(anno())[1+i]
    }
    selectInput(session$ns("voi"), "Choose the variable of interest", choices = variable_anno)
  })
  

  return(list(min_nodes_network = reactive({input$min_nodes_network}), 
              min_nodes_cluster = reactive({input$min_nodes_cluster}), 
              top_var = reactive({input$top_var}), 
              voi = reactive({input$voi}), 
              min_corr = reactive({input$min_correlation}), 
              range_corr = reactive({input$range_correlation}), 
              range_gfc = reactive({input$range_gfc})))  
}
