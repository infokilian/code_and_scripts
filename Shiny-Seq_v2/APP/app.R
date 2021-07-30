source('library.R')
source("module1.R")
source("PCA.R")
source("Boxplot.R")
source("module2.R")
source("remove_samples.R")
source("QC_normalized_data.R")
source("module3.R")
source("functions.R", encoding="UTF-8")
source("Differential_expression_project.R")
source("MA_plot.R")
source("Volcano_plot.R")
source("WGCNA.R")
source("CoCena_settings.R")
source("CoCena_data_processing.R")
source("CoCena_data_processing2.R", encoding="UTF-8")
source("CoCena_network_generation.R")
source("Kegg_module.R")
source("Biological_Process_module.R")
source("Hallmark_module.R")
source("FC-FC_plot_module.R")
source("Venn_diagram_module.R")
source("Enriched_markers_module.R")
source("Heatmap_module.R")
source("ANOVA_module.R")
source("TF_ChEA3_module.R")
source("powerpoint_module.R")

jsResetCode <- "shinyjs.reset = function() {history.go(0)}" # Define the js method that resets the page

ui<-tagList(
    navbarPage(
          title="Shiny-Seq",
          position="fixed-top",
          theme="design.css",
          id="ShinySeq",
               navbarMenu("Raw Data",
                       tabPanel("Raw Data",
                                 Module_Raw_data_UI("module"),
                                br(),
                                 useShinyjs(),# Set up shinyjs
                                 # Add a CSS class for red text colour
                                 inlineCSS(list(.blue = "background: lightblue")),
                                  uiOutput("button"),
                                br(),
                                 #Alter user if the number of samples in count/expression data
                                 #is not same as annotation table
                                fluidRow(column(1,downloadButton('anno_data', 'Download annotation')),
                                         hr(),
                                      column(1,downloadButton('exp_data', 'Download count data'))),
                                
                                
                                 fluidRow(column(8,
                                                bsAlert("alert_app"))),
                                br(),
                                br(),
                                hr(),
                                 #Display annotation table (excluded due to download from table directly)
                                 # downloadButton('downloadepData', 'Download full Data'),
                                 fluidRow(column (10,DT::dataTableOutput("pData"))),
                                br(),
                                br(),
                                hr(),
                                 #Display expression table
                                 fluidRow(
                                   column(1,
                                          selectInput(inputId= "datachoice"  ,label = h5("Select Data Type"), 
                                                      choices = list("Excel" = 1, "CSV" = 2),
                                                      selected = 1)),
                                   column(1, 
                                          br(),
                                          br(),
                                          downloadButton('downloadeData', 'Download Expression Data Table'))),
                                          
                                 fluidRow(column(12,DT::dataTableOutput("edata"))
                                 ),
                                 fluidRow(column(11, offset = 11,
                                                 br(),
                                                 textOutput("counter")))),
                               

                        tabPanel("Boxplot",
                                 Boxplot_module_UI("module")
                          ),
                       tabPanel( "PCA",
                                 br(),
                                 br(),
                                 br(),
                                 PCA_UI("module")
                       )
               ),
               navbarMenu("Normalization",
                          tabPanel("Normalized table",
                                   br(),
                                   br(),
                                    conditionalPanel(condition="input.ok1==0",                          #                           
                                                     p("please press start in the unnormalized table tab")),
                                    conditionalPanel(condition="input.ok1 >0",
                                                     fluidRow(column(4,uiOutput("condition"))),                 
                                                     fluidRow(column(5, uiOutput("design"))),
                                                     br(),
                                                     useShinyjs(), # Set up shinyjs
                                                    # Add a CSS class for red text colour
                                                     inlineCSS(list(.blue = "background: lightblue")),
                                                     uiOutput("button1"),
                                                    br(),
                                                    br(),
                                                    actionButton("rem_samp", "Detect outliers"),
                                                    actionButton("qc","Quality control"),
                                                    bsModal("modalqc","QC for normalized data","qc",size = "large",
                                                            conditionalPanel("input.ok2==0",
                                                                             p("Please press start button in the normalized tab and wait for a few minutes prior to clicking this button")),
                                                            conditionalPanel("input.ok2>0",
                                                                             QC_normalized_UI("module"))),
                                                    Module_Normalized_data_UI("module"),
                                                    
                                                    bsModal("modaloutlier", "Detect and remove outlier samples in addition to identifying batch effects", "rem_samp",size = "large",#uiOutput("plots")),
                                                            conditionalPanel("input.ok2==0",
                                                                             p("please press start button in the normalized tables tab and wait for a few minutes
                                                                               prior to clicking this button")
                                                                             ),
                                                            conditionalPanel("input.ok2>0",
                                                            Remove_samples_UI("module"),
                                                            bscols(column(3,actionButton("del", "Delete", icon = icon("trash-alt")),actionButton("reset", "Reset"))),
                                                            bscols(
                                                              column(10,DT::dataTableOutput("outliers")))
                                                            
                                                    ))

                         )
               ),
                          tabPanel("Boxplot",
                                   conditionalPanel(condition = "input.ok2==0",
                                                    p("Please press start button in normalized table tab")),
                                   conditionalPanel(condition="input.ok2 >0",
                                   Boxplot_module_UI("module1"))
                          ),
               tabPanel("Exploratory data anlysis",
                        tabsetPanel(
                           tabPanel("PCA",
                                    conditionalPanel(condition = "input.ok2==0",
                                                     p("Please press start button in normalized table tab")),
                                    conditionalPanel(condition = "input.ok2>0",
                                                     PCA_UI("moduleX"))),
                            tabPanel("Heatmap of samples",
                                     conditionalPanel(condition = "input.ok2==0",
                                                      p("Please press start button in normalized table tab")),
                                     conditionalPanel(condition = "input.ok2>0",
                                                     plotOutput("heatmap_sample")))
                        )
               )  
               ),
          tabPanel("Batch effect analysis",
                   tabsetPanel(
                     tabPanel("Compute batch effect",
                              br(),
                              actionButton("ok3","Start", icon = icon("paper-plane")),
                              conditionalPanel("input.ok3>0",
                                               Module_Batch_effect_UI("module"))
                              ),
                     
                     tabPanel("Source of variation",
                              conditionalPanel(condition = "input.ok3==0",
                                               p("Please press the start button in the batch-correction tab")),
                              conditionalPanel(condition = "input.ok3>0",
                                               uiOutput("sov_data_input"),
                                               uiOutput("go_sov"),
                                               plotlyOutput("sov")
                                               )),
                     tabPanel("Heatmap of samples after batch correction ",
                                       conditionalPanel(condition = "input.ok3==0",
                                                        p("Please press the start button in the batch-correction tab")),
                                       conditionalPanel(condition = "input.ok3>0",
                                                        plotOutput("heatmap_sample_batch"))
                              ))),
         
          navbarMenu("Downstream Analysis",
                       tabPanel("Differential expression analysis",
                                tabsetPanel(
                                  tabPanel("Get list of DE genes",
                                           conditionalPanel("input.ok3>0",
                                                            Module_Differential_Expression_UI("module")),
                                           conditionalPanel("input.ok3==0",
                                                            p("Please press start button in batch effect analysis tab"))
                                           ),
                                  tabPanel("Get list of Transcription factors",
                                           Enriched_markers_module_UI("module")
                                           
                                  ),
                                  tabPanel("Get list of marker genes",
                                           helpText("This tab enables you to search for any type of markers
                                                                 in the file that are enriched in the dataset uploaded.
                                                    Please make sure that the file only contains one column"),
                                           fluidRow(
                                             column(8,
                                                    fileInput("add_file", label = h3('Choose file to upload'),
                                                              accept = c(
                                                                'text/csv',
                                                                'text/comma-separated-values',
                                                                'text/tab-separated-values',
                                                                'text/plain',
                                                                '.csv',
                                                                '.tsv')
                                                    ))),
                                           Enriched_markers_module_UI("module1")
                                           
                                  ),
                                  tabPanel("Transcription factor prediction",
                                           predicted_TF_module_UI("module")
                                  ),
                                  tabPanel("Visualization",
                                           tabsetPanel(
                                              tabPanel("MA plot",
                                                       MA_plot_UI("module")
                                                       ),
                                              tabPanel("Volcano plot",
                                                       Volcano_plot_UI("module")),
                                              tabPanel("Venn Diagram",
                                                       Venn_diagram_module_UI("module")),
                                              tabPanel("FC-FC plot",
                                                       FC_FC_plot_module_UI("module"))
                                           ))
                                
                                )
                       ),
                     tabPanel("ANOVA",
                              tabsetPanel(
                                tabPanel("ANOVA Table",
                                         ANOVA_module_UI("module")
                                         ),
                                tabPanel("heatmap of 1000 most variable genes",
                                         heatmap_module_UI("module3")
                                )
                              )
                              
                     ),
                     tabPanel("Heatmap",
                              tabsetPanel(
                                tabPanel( "Heatmap of DE genes",
                                                           heatmap_module_UI("module")
                                          ),
                                tabPanel("heatmap of TF",
                                         helpText("This module is only available for the organisms 'Homo sapiens' and 'Mus musculus'"),
                                                          heatmap_module_UI("module1")
                                         ),
                                tabPanel("heatmap of marker genes",
                                         heatmap_module_UI("module2")
                                         )
                                
                                )
                              ),
                     tabPanel("Enrichment Analysis",
                              tabsetPanel(
                                tabPanel("Enrichment analysis: KEGG",
                                         Kegg_module_UI("module")),
                                tabPanel("Enrichment analysis-GO:Biological processes",
                                         Biological_process_module_UI("module")),
                                tabPanel("Hallmark plot",
                                         Hallmark_module_UI("module"))
                              ))
                     ),
          tabPanel("CoCena",
                   tabsetPanel(
                     tabPanel("1. Network settings",
                              conditionalPanel("input.ok3>0",
                                               CoCena_module_UI("module")),
                              br(),
                              conditionalPanel("input.ok3>0",
                                               actionButton("CoCenaButton","Save global settings", icon = icon("save")))),
                     tabPanel("2. Data processing I",
                              conditionalPanel("input.CoCenaButton>0",
                                               CoCena_data_processing_UI("module")),
                              conditionalPanel("input.CoCenaButton>0",
                              fluidRow(column(1,actionButton("confirm_choice", "Confirm cut-off", icon = icon("check-circle"))))),
                              br()),
                     tabPanel("3. Data processing II",
                              conditionalPanel("input.confirm_choice>0",
                                               CoCena_data_processing2_UI("module")),
                              fluidRow(column(1,actionButton("calc_network", "Calculate network", icon = icon("paper-plane")),offset=10))),
                     
                     tabPanel("4. Network generation",
                              conditionalPanel("input.calc_network>0",
                                               CoCena_network_generation_UI("module")))
                   )),
          
          tabPanel("Summary and results ", icon = icon("file-download"),
                   powerpoint_module_UI("module")
          )
          
                   
  )

)
#################################Begin of code#################
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 19MB.
options(shiny.maxRequestSize=500*1024^2)

server<-function(input, output,session) {
  
  output$counter <- 
    renderText({
      if (!file.exists("counter.Rdata")) counter <- 0
      else load(file="counter.Rdata")
      counter  <- counter + 1
      save(counter, file="counter.Rdata")     
      paste0("Visitors: ", counter)
    })
  
  #Delete all plots generated for ppt in the previous session
  do.call(file.remove,list(list.files("./plots",full.names = T)))
  input_data<-callModule(Module_Raw_data_Input,"module")
  print("inside app- line 98 module 1 pass")
  
  #display output (count data and annotation file)
  # 
  #V contains two reactive values wherein the variable v$data represents the count/expression table
  # v$data is updated when either of the buttons "delete" or "reset" is clicked
   v <- reactiveValues(data=NULL,click1=NULL,raw_pca=NULL,norm_pca=NULL,heat_de=NULL,heat_tf=NULL,heat_anova=NULL,
                       dds=NULL,batch=NULL,anova=NULL,kegg=NULL,bp=NULL,hallmark=NULL,
                       de=NULL,TF=NULL,overrep_TF=NULL,marker=NULL,data_wgcna=NULL,ma=NULL,vol=NULL,CoCena=NULL,
                       pheno_wgcna=NULL,
                       wgcna_output=NULL,
                       wgcna_click=FALSE)
   #makeReactiveBinding("v")
   
   ex_anno <- read.table("./example_data/annotation file.txt", header = TRUE,sep = "\t",check.names = FALSE)
   ex_exp <- read.table("./example_data/expression data.txt", header = TRUE,sep = "\t",check.names = FALSE)
   
   output$anno_data <- downloadHandler(
     
     filename = function(){
       paste("annotation file.txt")
     },
     content = function(file) {
       write.table(ex_anno, file = file,quote = FALSE, row.names = FALSE, sep = "\t")
     })  
   
   output$exp_data <- downloadHandler(
     
     filename = function(){
       paste("expression data.txt")
     },
     content = function(file) {
       write.table(ex_exp, file = file,quote = FALSE, row.names = FALSE, sep = "\t")
     })  


observeEvent(input$del, {
  print("inside app line 149")
  idx<-input$outliers_rows_selected
  data<-input_data$edata()
  data2 <-input_data$edata2()
  pheno<-input_data$pData()
  #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]
  data=data[,as.vector(pheno[,1])]
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Processing data", value = 0)
  # Increment the progress bar, and update the detail text.
  progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)

  if(!is.null(idx))
  {
    if(is.null(v$click1))
    {
      print("inside app line 181, delete button clicked hence deleteing selected samples")

      if(!identical(data, data2)) {
      data2$counts <- data2$counts[,-idx]
       data2$abundance <- data2$abundance[,-idx]
       data2$length <- data2$length[,-idx]}
      else{
        data2 <- data2[,-idx]
      }
      
      v$data<-list(data[,-idx],pheno[-idx,],rlogTransformation(as.matrix(data[,-idx])),data2)
      v$click1<-1
    }
    else
    {
      v$data[[1]]<-v$data[[1]][,-idx]
      v$data[[2]]<-v$data[[2]][-idx,]
      v$data[[3]]<-rlogTransformation(as.matrix(v$data[[1]][,-idx]))

    }
    v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1),
                      reactive({input_data$condition_module1()}),reactive({input_data$design_module1()}),reactive({input_data$filechoice()}))
    
    normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
    callModule(Boxplot_module,"module1",normal)
    callModule(PCA_module2,"moduleX",data=normal)
    v$norm_pca<-callModule(PCA_module2,"moduleX",data=normal)
    callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
    v_list<-callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                       v$dds$dds.fc()[[2]],
                       v$dds$normal())
    callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
    
  output$sov_data_input<-renderUI({
      checklist = list()
      for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
        
        checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
      
      selectizeInput("sov_inp", "Choose the variables to be included
                     to compute the amount of contribution to variance",
                     checklist,selected = NULL,multiple=T)
      
    })
  
  
  
    observeEvent(input$sov_inp,
                 {
                   output$go_sov<-renderUI({
                     actionButton("go_sov","Go")
                   })
                   observeEvent(input$go_sov,
                                {
                                  output$sov<-renderPlotly({
                                    id<-as.numeric(input$sov_inp)
                                    print("inside app line 528")
                                    if(length(id)!=0){
                                    df<-source_of_variation_op(v$dds$normal(),
                                                               colData(v$dds$dds.fc()[[1]]),
                                                               id)
                                    if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                    p <- plot_ly(df, labels = ~Variation,
                                                 values = ~Percentage, type = 'pie',
                                                 textposition = 'inside',
                                                 textinfo = 'label+percent',
                                                 insidetextfont = list(color = '#FFFFFF')) %>%
                                      layout(title = 'Source of variation plot',
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                    }
                                    else plotly_empty()
                                  })
                                })
                 })
      
    
    # Increment the progress bar, and update the detail text.
    progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
  }
})
#when reset button is clicked

observeEvent(input$reset, {
  #v$data <- rnorm(100)
  data<-input_data$edata()
  data2 <- input_data$edata2()
  pheno<-input_data$pData()
  #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Processing data", value = 0)
  # Increment the progress bar, and update the detail text.
  progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  v$data<-list(data,pheno,rlogTransformation(as.matrix(data)),data2)
  v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1),
                    reactive({input_data$condition_module1()}),reactive({input_data$design_module1()}),reactive({input_data$filechoice()}))
  
  normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
  callModule(Boxplot_module,"module1",normal)
  v$norm_pca<-callModule(PCA_module2,"moduleX",data=normal)
  callModule(PCA_module2,"moduleX",data=normal)
  callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
  callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                     v$dds$dds.fc()[[2]],
                     v$dds$normal())
  callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
  output$sov_data_input<-renderUI({

    checklist = list()
    for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
      
      checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
    
    selectizeInput("sov_inp", "Choose the variables to be included
                   to compute the amount of contribution to variance",
                   checklist,selected = NULL,multiple=T)
    
  })
  observeEvent(input$sov_inp,
               {
                 output$go_sov<-renderUI({
                   actionButton("go_sov","Go")
                 })
                 observeEvent(input$go_sov,
                              {
                                output$sov<-renderPlotly({
                                  id<-as.numeric(input$sov_inp)
                                  print("inside app line 499")

                                  if(length(id)!=0){
                                    edata<-v$dds$normal()
                                    if(v$batch$batch_choice()>1)
                                    {
                                      print("line 505 inside app")
                                      edata<-2^v$batch$batch_corrected_data()

                                    }
                                    df<-source_of_variation_op(v$dds$normal(),
                                                               colData(v$dds$dds.fc()[[1]]),
                                                               id)
                                    if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                    p <- plot_ly(df, labels = ~Variation,
                                                 values = ~Percentage, type = 'pie',
                                                 textposition = 'inside',
                                                 textinfo = 'label+percent',
                                                 insidetextfont = list(color = '#FFFFFF')) %>%
                                      layout(title = 'Source of variation plot',
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                    
                                  }
                                  else plotly_empty()
                                  
                                })
                              })
               })
  # Increment the progress bar, and update the detail text.
  progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
})  

  #Display expression table
  output$edata <- DT::renderDataTable({

    DT::datatable(v$data[[1]],class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 500,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list()
                  ))
    
  })
  
# download of expression table, either as excel or csv

  output$downloadeData <- downloadHandler(
    filename = function(){
      if(as.numeric(input$datachoice==1)){paste("Expression Table.xlsx")}
      else if (as.numeric(input$datachoice==2)){paste("Expression Table.csv")}
    },
    content = function (file){
      if(as.numeric(input$datachoice==1)){
  
        wb <- createWorkbook()
        addWorksheet(wb, sheetName = "Expression data")
        writeData(wb = wb, sheet = 1, x = v$data[[1]], colNames = T, rowNames = T)
        saveWorkbook(wb, file)

        }
      else if (as.numeric(input$datachoice==2)){write.csv(v$data[[1]], file)}
      
    }
  )   

  
   # #Display annotation table
  output$pData <- DT::renderDataTable(server=FALSE,{

    DT::datatable(v$data[[2]],class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download annotation data')#I('colvis')
                  ))
    )

  })


  #start button in raw data tab
  observeEvent(input_data$file2(),{
    if (!(as.numeric(input_data$refgenome()) == 1 & input_data$organism() %in% c("Anopheles gambiae",
                                                               "Arabidopsis thaliana",
                                                               "Canis familiaris",
                                                               "K-12",
                                                               "Sakai",
                                                               "Plasmodium falciparum",
                                                               "Others")))
    {
      closeAlert(session, "Use own upload")
      output$button<-renderUI({
        fluidRow(column(1,
        actionButton("ok1", label = "Start", icon("paper-plane"))),
        br()) 
      })
      } else {
        createAlert(session, "alert_app", "Use own upload", title = "Error",
                    content = "The chosen organism has no Ensembl repository. Please choose 'Own upload'")
                  
      }

  })

  #start button turns blue on click
  #call boxplot and PCA for raw data input
  observeEvent(input$ok1,{
    
    observeEvent(input_data$file2(),{
      
      if(!is.null(input_data$pData())) {
        closeAlert(session,"exampleAlert_app")
        print("inside app line 111")
        data<-input_data$edata()
        data2<-input_data$edata2()
        pheno<-input_data$pData()
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Processing Data", value = 0)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)

        data=data[,as.vector(pheno[,1])]
        print("line 291 app")
        #progress$set(message = "Rlog transformation step", value = 0)
        v$data<-list(data,pheno,as.matrix(data), data2)
        # Increment the progress bar, and update the detail text.
        progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      else {
        createAlert(session,"alert_app", "exampleAlert_app", title = "Oops!",
                    content = paste0("Some samples in the expression/count table is absent in the annotation table!
                  Please give correct input"), append = FALSE)
      }

    })
    
    toggleClass("ok1", "blue")
    callModule(Boxplot_module,"module",v$data)
    v$raw_pca<-callModule(PCA_module,"module",data=v$data) # save for ppt module
      callModule(PCA_module,"module",data=v$data) # show plot in app
 
    # creates a checkbox widget to select the important condition
    output$condition <-
      renderUI({
        if(as.numeric(input_data$filechoice())!=3){
        checklist = list()
        for (i in seq_along(colnames(v$data[[2]]))[-1]) {
          checklist[[colnames(v$data[[2]])[[i]]]] = i
        }
        radioButtons("conchoice", "Choose the conditon", checklist,selected = 2)
        }
      })
    
    # creates a checkbox widget to select the variables to be included in design
    output$design <-
      renderUI({
        if(as.numeric(input_data$filechoice())!=3){
        req(input$conchoice)
        checklist = list()
        for (i in seq_along(colnames(v$data[[2]]))[-1]) {
          # disply all variables except condition variable
          if(i!=as.numeric(input$conchoice)) checklist[[colnames(v$data[[2]])[[i]]]] = i
        }
        if(length(checklist)!=0){
          checkboxGroupInput(session$ns("designchoice"), "Choose the variables to be included in design", checklist)
        }
        }
      })
    output$button1<-renderUI({
      if(!is.null(input$conchoice)|!is.null(input_data$condition_module1())) actionButton("ok2", label = "Start", icon = icon("paper-plane"))
    })
    
  })

observeEvent(input$ok2,
             {
               toggleClass("ok2", "blue")
               print("inside app line 324 calling normalized table module")
               
               v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1),
                                 reactive({input_data$condition_module1()}),reactive({input_data$design_module1()}),reactive({input_data$filechoice()}))

               normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
               #normalized table,annotation table,rlog transformed data for PCA
               callModule(Boxplot_module2,"module1",normal)
               v$norm_pca<-callModule(PCA_module2,"moduleX",data=normal) # save for ppt module
                 callModule(PCA_module2,"moduleX",data=normal) # show plot in app
               callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)

               output$heatmap_sample<-renderPlot({
                 sampleDists <- dist(t(assay(v$dds$dds.fc()[[2]])))
                 print("inside app line 368")
                 heat_ok2 <- heatmap_sample(sampleDists,v$dds$dds.fc()[[2]])
                 heat_ok2
               })
               output$sov_data_input<-renderUI({
                     checklist = list()
                     for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
                       
                         checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
                     
                     selectizeInput("sov_inp", "Choose the variables to be included
                                      to compute the amount of contribution to variance",
                                    checklist,selected = NULL,multiple=T)
                   
                 })
               observeEvent(input$sov_inp,
                            {
                              output$go_sov<-renderUI({
                                actionButton("go_sov","Go")
                              })
                              observeEvent(input$go_sov,
                                           {
                                             output$sov<-renderPlotly({
                                               id<-as.numeric(input$sov_inp)
                                               print("inside app line 528")
                                               if(length(id)!=0){
                                                 edata<-v$dds$normal()
                                                 if(v$batch$batch_choice()>1)
                                                   {
                                                     edata<-2^v$batch$batch_corrected_data()
                                             
                                                 }
                                               df<-source_of_variation_op(edata,
                                                                          colData(v$dds$dds.fc()[[1]]),
                                                                                  id)
                                               if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                               p <- plot_ly(df, labels = ~Variation,
                                                            values = ~Percentage, type = 'pie',
                                                            textposition = 'inside',
                                                            textinfo = 'label+percent',
                                                            insidetextfont = list(color = '#FFFFFF')) %>%
                                                 layout(title = 'Source of variation plot',
                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                               }
                                               else plotly_empty()
                                             })
                                           })  
                            })
               
               
observeEvent(input$rem_samp,
             {
               v_list<-callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                                                          v$dds$dds.fc()[[2]],
                                                          v$dds$normal())
                                                          
               print("inside app line 328 bar")
               #display annotation table after removing outliers
               output$outliers <- DT::renderDataTable({
                 dat <- as.data.frame(colData(v$dds$dds.fc()[[2]]))
                 d <- v_list$d()
                 v1_clickBar_list<-v_list$bar()
                 idx<-NULL
                 if (!is.null(d) && is.null(v1_clickBar_list))
                 {
                   
                   idx<-which(dat[,1] %in% d)
                 }
                 else if(is.null(d) && !is.null(v1_clickBar_list))
                 {
                   idx<-round(v1_clickBar_list)
                 }
                 else if(!is.null(d) && !is.null(v1_clickBar_list))
                 {
                   id<-c(which(dat[,1] %in% d),round(v1_clickBar_list))
                   idx<-unique(id)
                 }
                 
                 library(DT)
                 DT::datatable(dat,class = 'cell-border stripe',
                               selection = list(mode='multiple',selected=idx,target='row'),
                               extensions = list('Scroller'=NULL,'Buttons'=NULL),
                               options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 350,scroller = TRUE,dom = 'Bfrtip',
                                              buttons = list()))
              })

             })

callModule(powerpoint_module,"module",
           reactive({v$dds$normal()}),reactive({v$dds$dds.fc()[[1]]}),reactive({input$conchoice}),
           reactive({v$raw_pca()}),reactive({v$norm_pca()}),reactive({normal}),reactive(v$data),
           reactive({v$batch$top_batch_pca()}),reactive({v$batch$batch_corrected_data()}),reactive({v$batch$batch_choice()}),
           reactive({v$de$combination()}),reactive({input$ok2}),reactive({input$ok3}),reactive({v$de$ok3()}),reactive({v$de$de_genes()}),reactive({v$TF$de_genes()}),
           reactive({v$ma$input_scale()}),reactive({v$ma$p_values()}),
           reactive({v$vol$input_scale_volx()}),reactive({v$vol$input_scale_voly()}),reactive({v$vol$hypothesis_choice()}),
           reactive({v$heat_de$input_Distance()}),reactive({v$heat_de$input_Linkage()}),
           reactive({v$heat_tf$input_Distance()}),reactive({v$heat_tf$input_Linkage()}),
           reactive({v$anova$dds()}),reactive({v$heat_anova$input_Distance()}),reactive({v$heat_anova$input_Linkage()}),
           reactive({v$kegg$Enriched_Kegg_table()}),reactive({v$kegg$Enriched_Kegg_obj()}),
           reactive({v$bp$Enriched_BP_table()}),reactive({v$bp$Enriched_BP_obj()}),
           reactive({v$hallmark$Enriched_hall_table()}),reactive({v$hallmark$Enriched_hall_obj()}),
           reactive({v$anova$anova_table()}), reactive({input_data$organism()}),
           reactive({input_data$dataset()}),reactive({v$batch$batch_data_for_DESeq()}),
           reactive({input_data$condition_module1()}),reactive({input_data$filechoice()})
)

observeEvent(input$ok3,
             {

               print("inside app line 406")
               v$batch<-callModule(Module_Batch_effect,"module",reactive({v$dds$conchoice_app()}),input$designchoice,reactive({v$dds$dds.fc()}),reactive({v$dds$normal()}))

               normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))


               v$CoCena$settings<-callModule(CoCena_module,"module",
                                    reactive({input_data$pData()}))
               
               observeEvent(input$CoCenaButton,{
               v$CoCena$data_processing1<-callModule(CoCena_data_processing,"module",
                                                    reactive({v$CoCena$settings$top_var()}),
                                                    reactive({v$CoCena$settings$range_corr()}),
                                                    reactive({v$CoCena$settings$min_corr()}),
                                                    reactive({v$CoCena$settings$min_nodes_network()}),
                                                    reactive({v$batch$batch_choice()}),
                                                    reactive({v$batch$batch_corrected_data()}),
                                                    reactive({v$dds$normal()})
               )
               
               v$CoCena$data_processing2<-callModule(CoCena_data_processing2,"module",
                                                     reactive({input_data$pData()}),
                                                     reactive({v$CoCena$data_processing1$top_var_output()}),
                                                     reactive({v$CoCena$data_processing1$cutoff_choice()}),
                                                     reactive({v$CoCena$data_processing1$gene_corr()}),
                                                     reactive({v$CoCena$settings$voi()}),
                                                     reactive({v$CoCena$settings$range_gfc()})
                                                     )
               v$CoCena$network<-callModule(CoCena_network_generation,"module",
                                            reactive({v$CoCena$data_processing2$filt_cutoff()}),
                                            reactive({v$CoCena$data_processing2$GFC_all_genes()}),
                                            reactive({v$CoCena$settings$min_nodes_cluster()}),
                                            reactive({v$CoCena$settings$voi()})
                                            )
               })
              
               v$de<-callModule(Module_Differential_Expression,"module",input$conchoice,
                                reactive({v$batch$batch_data_for_DESeq()}),
                                reactive({v$CoCena$network}),
                                reactive({v$dds$normal()}),
                                reactive({v$batch$batch_choice()}),
                                reactive({v$batch$batch_corrected_data()}),
                                reactive({v$anova$anova_table()}),
                                reactive({input_data$condition_module1()}),
                                reactive({input_data$filechoice()})
                                )

               callModule(FC_FC_plot_module,"module",reactive({v$de$de_genes}),
                          reactive({v$de$combination})
                          )

               v$ma<-callModule(MA_plot,"module",reactive({v$de$combination}),
                          reactive({v$de$de_genes()}),
                          reactive({v$de$p_value()}))
               v$vol<-callModule(Volcano_plot,"module",
                          reactive({v$de$combination()}),
                          reactive({v$de$de_genes()}),
                          reactive({v$de$p_value()}),
                          reactive({v$de$hypothesis_choice()}),
                          reactive({input_data$organism()}),
                          reactive({input_data$dataset()}))
               
               callModule(Venn_diagram_module,"module",reactive({v$de$de_genes()}),
                                     reactive({v$de$combination}),
                                     reactive({v$CoCena$network}))
               
                             v$heat_de<-callModule(heatmap_module,"module","DE",NULL,
                                        reactive({v$batch}),
                                        reactive({v$de$de_genes()}),
                                        reactive({v$de$combination()}),
                                        reactive({v$CoCena$network}),
                                        reactive({input_data$organism()})

                             )
 
               observeEvent(input$add_file,
                            {
                              filepath<-input$add_file$datapath
                              TF_list<-read.csv(filepath, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                              
                              v$marker<-callModule(Enriched_markers_module,"module1","marker", reactive({TF_list}),
                                                 reactive({input_data$organism()}),
                                                 reactive({v$de$de_genes()}),
                                                 reactive({v$de$combination}),
                                                 reactive({v$CoCena$network}),
                                                 reactive({v$dds$conchoice_app()}),
                                                 reactive({v$dds$normal()}),
                                                 reactive({v$batch$batch_data_for_DESeq()}),
                                                 reactive({v$batch$batch_choice()}),
                                                 reactive({v$batch$batch_corrected_data()}),
                                                 reactive({v$anova$anova_table()}),
                                                 reactive({v$dds$dds.fc()[[1]]})
                                                 )
                              
                          v$heat_marker<-callModule(heatmap_module,"module2","TF_custom",NULL,
                                         reactive({v$batch}),
                                         reactive({v$marker$de_genes()}),
                                         reactive({v$de$combination()}),
                                         reactive({v$CoCena$network}),
                                         reactive({input_data$organism()}) 
                              )
                            })

                                     
               TF_list<-read.csv("./www/Transcriptome_TFcat.txt", header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
               if (input_data$organism()=="Homo sapiens"){ 
                 TF_list <- TF_list[,"Human"]
                 
               } else if (input_data$organism()=="Mus musculus"){ 
                 TF_list <- TF_list[,"Mouse"]
                 
               }
          v$TF<-callModule(Enriched_markers_module,"module",
                           "TF",
                          reactive({TF_list}),
                          reactive({input_data$organism()}),
                          reactive({v$de$de_genes()}),
                          reactive({v$de$combination}),
                          reactive({v$CoCena$network}),
                          reactive({v$dds$conchoice_app()}),
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_data_for_DESeq()}),
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_corrected_data()}),
                          reactive({v$anova$anova_table()}),
                          reactive({v$dds$dds.fc()[[1]]})
                          )
          
          v$heat_tf<-callModule(heatmap_module,"module1","TF",NULL,
                     reactive({v$batch}),
                     reactive({v$TF$de_genes()}),
                     reactive({v$de$combination()}),
                     reactive({v$CoCena$network}),
                     reactive({input_data$organism()})
                       
          )
               
           v$anova<-callModule(ANOVA_module,"module",
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$batch$batch_data_for_DESeq()[[2]]}),
                          reactive({v$de$combination()}),
                          reactive({v$dds$conchoice_app()}),
                          reactive({v$de$de_genes}),
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_corrected_data()}),
                          reactive({v$de$ok3()})
                
               )
           
           
           v$heat_anova<-callModule(heatmap_module,"module3","ANOVA",
                                 reactive({v$anova$dds()}),
                                 reactive({v$batch}),
                                 reactive({v$de$de_genes()}),
                                 reactive({v$de$combination()}),
                                 reactive({v$CoCena$network}),
                                 reactive({input_data$organism()})
                                 
           )

            callModule(predicted_TF_module,"module",
                          reactive({v$de$de_genes()}),
                          reactive({v$batch$batch_data_for_DESeq()}),
                          reactive({v$anova$anova_table()}),
                          reactive({v$de$combination()}),
                          reactive({v$dds$conchoice_app()}),
                          reactive({v$CoCena$network}),
                          reactive({input_data$organism()}),
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_corrected_data()}),
                          reactive({input_data$dataset()})
                          
               )
               v$kegg<-callModule(Kegg_module,"module",reactive({v$de$de_genes}),
                                               reactive({input_data$organism()}),
                                               reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                                               reactive({v$de$combination}),
                                               reactive({v$CoCena$network}),
                                               reactive({v$anova$anova_table()}),
                                               reactive({input_data$orgDB()})
                          )
                              
               v$bp<-callModule(Biological_process_module,"module",reactive({v$de$de_genes}),
                          reactive({input_data$organism()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$de$combination}),
                          reactive({v$CoCena$network}),
                          reactive({input_data$orgDB()}))
               
               v$hallmark<-callModule(Hallmark_module,"module",reactive({v$de$de_genes}),
                          reactive({input_data$organism()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$de$combination}),
                          reactive({v$CoCena$network}),
                          reactive({input_data$orgDB()}))
               
               print("inside app line 455")
              
               output$heatmap_sample_batch<-renderPlot({
                 sampleDists <- dist(t(assay(v$batch$batch_data_for_DESeq()[[3]])))
                 var<-as.numeric(v$dds$conchoice_app())
                 heat_batch <- heatmap_sample(sampleDists,v$batch$batch_data_for_DESeq()[[3]])
                 heat_batch
               })
             })
            })

    
output$file_input1<-renderUI({
        req(input$option_wgcna)
  if(as.numeric(input$option_wgcna)==1)
  {
        fileInput("input_file",label=h5("Please upload normalized/batch corrected data"),
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'))
  }
})
        
  observeEvent(input$input_file,
                     {
                       filepath<-input$input_file$datapath
                       if(!is.null(filepath))
                       {
                         input_exp<-read.csv(filepath, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                         data<-input_exp[,-1]
                         
                         #convert entries to integers
                         data=as.matrix(data)
                         storage.mode(data)="double"
                         data= data.frame(data)
                         rownames(data)=input_exp[,1]
                         v$data_wgcna<-data
                         
                         output$file_input2<-renderUI({
                           fileInput("input_file2",label=h5("Please upload annotation data"),
                                     accept = c('text/csv',
                                                'text/comma-separated-values',
                                                'text/tab-separated-values',
                                                'text/plain',
                                                '.csv',
                                                '.tsv'))
                         })
                         
                         observeEvent(input$input_file2,
                                      {
                                        filepath2<-input$input_file2$datapath
                                        if(!is.null(filepath2))
                                        {
                                          input_p<-read.csv(filepath2, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                                          pheno<-input_p
                                          v$pheno_wgcna<-get_pheno(v$data_wgcna,pheno)
                                          
                                        }
                                        
                                      })
                       }
                       
                     })
      
  
output$pheno_table <- DT::renderDataTable({
  DT::datatable(v$pheno_wgcna,class = 'cell-border stripe',
                selection = list(mode='single',target = 'column'),
                extensions = list('Scroller'=NULL,'Buttons'=NULL),
                options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                               buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download')
                ))
  )
  
})
output$go_wgcna<-renderUI({
  req(input$option_wgcna)
  if(as.numeric(input$option_wgcna)==1)
  {
    req(input$pheno_table_columns_selected)
  if (as.numeric(input$pheno_table_columns_selected)>1) {actionButton("go_wgcna","GO")} #any variable other than sample id(assuming saple id to be the first column)
  }
  else if(as.numeric(input$option_wgcna)==2) actionButton("go_wgcna","GO")
  
})
observeEvent(input$pheno_table_columns_selected,{
  col<-as.numeric(input$pheno_table_columns_selected)
  condition<-v$pheno_wgcna[,col]

  
  observeEvent(input$go_wgcna,
               {
                 callModule(WGCNA_module,"module",
                            infile=reactive({v$data_wgcna}),perform_voom=FALSE,condition,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            reactive({0}),
                            NULL,
                            NULL
                            )
                 
               })
  
})

}
shinyApp(ui = ui, server = server)

