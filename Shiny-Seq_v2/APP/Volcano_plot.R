source("functions.R")
Volcano_plot_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    fluidRow(column(3,uiOutput(ns("vol_comb"))),
             column(3,textInput(ns("max_number"), "Choose how many genes should be annoted:", value = 10))),
    br(),
    fluidRow(column(3,uiOutput(ns("vol_y"))),
             column(3,uiOutput(ns("anno_hallmark"))),
             column(3,radioButtons(ns("user_choice"), "Custom annotation", choices=c("Specified genes", "Upload gene list"), selected = "Specified genes"))),
    br(),
    fluidRow(column(3, uiOutput(ns("vol_x"))),
             column(3,uiOutput(ns("anno_GO"))),
             column(3,uiOutput(ns("custom")))),
    br(),
    fluidRow(column(2, actionButton(ns("Volcano"),"Go!"))),
    br(),
    plotOutput(ns("volcano_plot"), width = "800px", height = "600px")
    
  )
}
Volcano_plot<-function(input,output,session,combination,DE_genes,p_values,hypothesis_choice,organism, dataset)
{
  
  output$custom <- renderUI({
    if(input$user_choice == "Specified genes"){
      textInput(session$ns("custom_gene"), "Annotate specified genes", placeholder = "i.e. A2M, ABL1, ADCY5, etc.", value ="")
    } else if(input$user_choice == "Upload gene list"){
      
      fileInput(session$ns("custom_list"), "Annotate uploaded gene list",
                accept = c('.csv',
                           '.tsv',
                           'text/comma-separated-values',
                           'text/csv',
                           'text/plain',
                           'text/tab-separated-values'))
    }
  })
  
  user_gene <- reactive({
    if(input$user_choice == "Specified genes") input$custom_gene
    else NULL
  })
  
  user_list <- reactive({
    if(input$user_choice == "Upload gene list") input$custom_list
    else NULL
  })

  #combinations for volcano plot (A vs B, C vs D)
  output$vol_comb <- renderUI({
    if(!is.null(combination()))
    {
      num<- length(combination())
      comb<-lapply(1:num, function(i) {
        input$combination[i]
      })
      
      checklist<-list()
      for (i in seq_along(comb)) {
        checklist[[comb[[i]]]] = i
      }
      selectInput(session$ns("vol_choice"),label = h5("Choose comparison") ,
                  choices = checklist,selected = 1)
    }
  })
  
  #Slider to adjust Y axis of volcano plot
  output$vol_y<-renderUI({
    if(!is.null(input$vol_choice))
    {
      sliderInput(session$ns("scale_voly"),"p-value cut-off",0,0.5,value = 0.05,step = 0.001)
    }
  })
  #slider to adjust x axis of volcano plot
  
  output$vol_x<-renderUI({
    if(!is.null(input$vol_choice))
    {
      num<- length(combination())
      #Get the deseq2 dataset
      res<-DE_genes()[[as.numeric(input$vol_choice)]][5] [[1]]
      max<-max(na.omit(abs(res$log2FoldChange)))
      max<-round(max,1)
      sliderInput(session$ns("scale_volx"),"log2(FoldChange) cut-off",0,max+5,value =1,step = 0.1)
    }
  })
  output$anno_hallmark<-renderUI({
    if (organism() %in% c("Homo sapiens", 
                          "Mus musculus", 
                          "Bos taurus", 
                          "Caenorhabditis elegans", 
                          "Canis familiaris",
                          "Danio rerio",
                          "Drosophila melanogaster",
                          "Gallus gallus",
                          "Rattus norvegicus",
                          "Sus scrofa")){
    if(!is.null(input$vol_choice))
    {
      if (organism() == "Canis familiaris")  c1_hallmark <- msigdbr(species = "Canis lupus familiaris", category = "H")
      else c1_hallmark <- msigdbr(species = organism(), category = "H")
      gene_sets <- list()
      for (p in 1:(length(unique(c1_hallmark))+1)){
        if(p==1){
          gene_sets[[p]] <- "Default"
          names(gene_sets)[[p]] <- "Default"
        }else{
        gene_sets[[p]] <- unique(c1_hallmark$gs_name)[p-1]
        names(gene_sets)[[p]] <- unique(c1_hallmark$gs_name)[p-1]
        }
      }
      selectInput(session$ns("hallmark"),"Annotate Hallmark gene set",choices = gene_sets)
    }
  }
  })
  
  output$anno_GO<-renderUI({
    if(!(organism() %in% c("Anopheles gambiae",
                         "Arabidopsis thaliana",
                         "Canis familiaris",
                         "K-12",
                         "Sakai",
                         "Plasmodium falciparum",
                         "Others"))){
    if(!is.null(input$vol_choice))
    {
   
      textInput(session$ns("GO"),"Annotate GO term", placeholder = "i.e. GO:0031981", value ="")
    }
    }
  })
  
  #volcano plot
  volcano_seq_plot <- eventReactive(input$Volcano,{
    print("inside volcano plot line 64")
    req(input$vol_choice)
    req(input$scale_volx)
    req(input$scale_voly)
    if(!is.null(input$vol_choice) & !is.null(input$scale_volx) & !is.null(input$scale_voly))
    {
      volcano_seq <-volcano_plot(input$vol_choice,combination(),DE_genes(),input$scale_volx,input$scale_voly,
                                 input$hallmark,input$GO,organism(),dataset(),user_gene(),user_list(),as.numeric(input$max_number))
      volcano_seq
    }
  }) 

  output$volcano_plot<- renderPlot({
   
      volcano_seq_plot()

  })
  

  
return(list(
  input_scale_volx=reactive({input$scale_volx}),
  input_scale_voly=reactive({input$scale_voly}))
)  
}