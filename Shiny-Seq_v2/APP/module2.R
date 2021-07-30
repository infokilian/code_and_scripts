source("functions.R")
source("gene_count_module.R")
Module_Normalized_data_UI<-function(id)
{
  ns<-NS(id)
  tagList(
          fluidRow(column(10,bsAlert("message"))),
          fluidRow(column(5,actionButton(ns("help"), "Help", icon = icon("question-circle")))),
          fluidRow(column(8,
                            bsAlert("alert"),
                            plotlyOutput(ns("cut_off_plot"))
                          ),
                   column(4,helpText("Please enter before pressing Start"),
                          textInput(ns("cutoff"), "Cut-off value", value = 10, width = NULL, placeholder = NULL))),
          fluidRow(
            column(1,
                   selectInput(ns("datachoice2")  ,label = h5("Select Data Type"), 
                               choices = list("Excel" = 1, "CSV" = 2),
                               selected = 1)),
            column(1, 
                   br(),
                   br(),
                   downloadButton(ns("downloadnorm"), 'Download normalized Data Table'))),
          fluidRow(
            column(12,
                 DT::dataTableOutput(ns("normal")),
                 gene_count_module_UI(ns("module"))
          )),
          
          bsModal(ns("modalhelp"), "Help for error model matrix not full rank", ns("help"),size = "large",
                             helpText("There are two main reasons for this problem: either one or more
                                      columns in the model matrix are linear combinations of other columns, or there
                                      are levels of factors or combinations of levels of multiple factors which are missing
                                      samples. We address these two problems below"),
                             helpText("Case 1:"),
                             fluidRow(column(5,img(src="error1.png"))),
                             helpText(" In the above table two variables contain exactly the same information. The software cannot
                                      fit an effect for batch and condition, because
                                      they produce identical columns in the model matrix. This is also referred to as
                                      perfect confounding."),
                             p("Solution:"),
                             helpText("the batch effect cannot be fit and must be removed
                                      from the model formula. There is just no way to tell apart the condition effects
                                      and the batch effects. The options are either to assume there is no batch effect
                                      (which we know is highly unlikely given the literature on batch effects in sequencing
                                      datasets) or to repeat the experiment and properly balance the conditions across
                                      batches."),
                             helpText("Case 2:"),
                             fluidRow(column(5,img(src="error2.png"))),
                             helpText("In the above table the variables are not identical,
                                      but one variable can be formed by the combination of other factor levels. In the
                                      following example, the effect of batch 2 vs 1 cannot be fit because it is identical
                                      to a column in the model matrix which represents the condition C vs A effect."),
                             p("Solution:"),
                             helpText("the batch effect cannot be fit and must be removed
                                      from the model formula. There is just no way to tell apart the condition effects
                                      and the batch effects. The options are either to assume there is no batch effect
                                      (which we know is highly unlikely given the literature on batch effects in sequencing
                                      datasets) or to repeat the experiment and properly balance the conditions across
                                      batches."),
                             helpText("Case 3:"),
                             fluidRow(column(5,img(src="error3.png"))),
                             helpText("In the above case where we can in fact perform inference.In the above example the experiment
                                      comprises of grouped individuals, where the goal is to test the group-specific effect of
                                      a treatment, while controlling for individual effects."),
                             p("Solution:"),
                             helpText("unfortunately the tool cannot be currently used for a similar scenario
                                      as above.We are working on it.You will have to use DESeq2 to define the design."),
                             helpText("Case 4:"),
                             fluidRow(column(5,img(src="error4.png"))),
                             helpText("In the above case as highlighted a level is missing from a factor.The group 3 does not have level C"),
                             p("Solution:"),
                             helpText("Remove these samples with missing level.In the above case remove samples 13 to 16")
                           )#,
               
    )
  
}
#normalized tab
#create a gene plots module
Module_Normalized_data<-function(input,output,session,vdata,conchoice,designchoice,v_click1,conchoice_module1,designchoice_module1,filechoice)
{
  data<-vdata()
  trim <- reactive({
    tryCatch(
      
      #check if user has provided the choice of treatment/condition variable
       {
        print("treatment/condition variable provided")
        #Get variable to consider as treatment/condition
        if (as.numeric(filechoice())==3) condition <- as.numeric(conchoice_module1()) #already chosen for STAR
        else condition<-as.numeric(conchoice)

        
        #===============Main purpose of this reactive==========================#
        #Step 1: prepare for DESeq2
        #Inputs needed: expression/count table, annotation table and design formula
        
        #Step 2: prepare expression/count table for trimming (Quality control used to remove genes with low counts)
        #=======================================================================#
        
        #Input 1: Get expression data
        edata <-data[[4]]
        
        #Input 2: Get annotation data
        pheno<-data[[2]]

        #Input 3: Get design formula
        if (as.numeric(filechoice())==3) design<-as.numeric(designchoice_module1()) #already chosen for STAR
        design<-as.numeric(designchoice)

        #Rename treatment/condition variable to condition
        #This is done to ensure one standard name to refer to , in the future computation steps
        colnames(pheno)[condition]<-"condition"
        
        #Since the design formula contains the name of the treatment/condition variable(based on input annotation table)
        #We rename the variable to "condition" in the design formula
        #Updating design formula
        d<-''
        for (i in design)
        {
          d<-paste(d,colnames(pheno)[i]," + ",sep=" ")
        }
        d<-paste('~',' ',d,colnames(pheno)[condition],sep = " ")

        #once we have obtained the input we store it in a container (a subclass)
        #DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values,
        #intermediate calculations and results of an analysis of differential expression.
        #check dese2 vignette for more information.

        #Prepare DESeq dataset using inputs

        if(identical(data[[4]], data[[1]]))
          dds <- DESeqDataSetFromMatrix(edata,colData=pheno,formula(d))
        
        else
          dds <- DESeqDataSetFromTximport(edata,colData=pheno,formula(d))
        
        #A check to see if the annotation table has been correctly stored and if the value can be retrieved

        genes_to_keep <- rowSums(counts(dds)) >= as.numeric(input$cutoff)
        dds <- dds[genes_to_keep,]
        dds <- DESeq(dds)
        print("dds")
        print(dds)
        print(colData(dds))
        
        dds.norm <- as.data.frame(counts(dds, normalized=T))
        
        if(identical(data[[4]], data[[1]]))
          dt1<-data.frame(rowSum = rowSums(edata))
        
        else
        dt1<-data.frame(rowSum = rowSums(edata$counts))
        #Normalization: DESeq2 performs normalization by first computing size factor for each sample
        
        #call function min_samples_three to get a list of treatment/condition groups that contain only one or two samples
        samp<-min_samples_three(dds.norm)
        
        #prior to proceeding with normalization, if treatment/condition groups have one or two samples
        #then user is alerted through a warning message (error message if a treatment/condition group has only one sample) 
        
        if(is.null(samp)) {
          closeAlert(session,"exampleAlert")
        }
        else if(samp[[2]]==1) { #get treatment/condition groups with only one sample
          if(length(samp[[1]])>1) #if there are more than one treatment/condition group with only one sample
          {
            x<-strsplit(samp[[1]],' ') #Get names of all treatment/condition groups with only one sample
            if(!is.null(x)) closeAlert(session,"exampleAlert")
            sam<-NULL
            for (i in 1:(length(x)-1))
            {
              sam<-paste0(sam,x[[i]],',')
            }
            sam<-paste0(substr(sam,1,nchar(sam)-1),' and ',x[[length(x)]])
            #Alert user with the  names of all treatment/condition groups with only one sample
            createAlert(session,"alert", "exampleAlert", title = "Oops!",
                        content = paste0("The conditions ",sam," have only 1 sample.Either add additional samples for these conditions or remove the condition."), append = T)
          }
          else{
            #Alert user with the  name of  treatment/condition group with only one sample
            createAlert(session,"alert", "exampleAlert", title = "Oops!",
                        content = paste0("The condition ",samp[[1]]," has only 1 sample.Either add additional samples for this condition or remove the condition."), append = T)
          }
          
          
        }
        else if(samp[[2]]==2) #get treatment/condition groups with only two samples
        {
          if(length(samp[[1]])>1)#if there are more than one treatment/condition groups with only two samples
          {
            x<-strsplit(samp[[1]],' ')#Get names of all treatment/condition groups with only two samples
            if(!is.null(v_click1())) print(x)
            if(!is.null(x)) closeAlert(session,"exampleAlert")
            sam<-NULL
            for (i in 1:(length(x)-1))
            {
              sam<-paste0(sam,x[[i]],',')
            }
            sam<-paste0(substr(sam,1,nchar(sam)-1),' and ',x[[length(x)]])
            
            #Alert users with names of all treatment/condition groups with only two samples
            createAlert(session,"alert", "exampleAlert", title = "Warning!",
                        content = paste0("The conditions ",sam," have only 2 samples"), append = T)
            
          }
          else
          {
            #Alert user with name of treatment/condition group with only two samples
            createAlert(session,"alert", "exampleAlert", title = "Warning!",
                        content = paste0("The condition ",samp[[1]]," has 2 samples"), append = T)
            
          }
          
        }
        list(dds, dds.norm, dt1)
  
    }
      #It is possible that two variable in the annotation table can coontain the same information
      #meaning the variables explain the same variance in the dataset in which case DESEq2 throws an error
      #"the model matrix is not full rank". In such a case the user may have to remove one of these
      #variables from the annotation table or combine it as one. For more information rgearding 
      #the problem and what to do please check for the error in DESeq2 vignette.
      , error = function(c) {
        print(typeof(c$message))
        print(substr(c$message,1,33))
        if("the model matrix is not full rank" == substr(c$message,1,33))  c$message <- paste0("Error:The model matrix is not full rank. To know more about this please click the help button on this page.")
        else c$message
        stop(c)
      })
    
  })
#Prior to trimming( remove genes with low counts), we need to get the threshold from the user
#The user decides the the threshold based on the cut off plot.
# The cut off plot is an interactive plot which indicates the number of genes remaining in the dataset after
# trimming. The user can hover over the plot to see the threshold value and the corresponding genes that would remain in the 
# dataset if the value is chosen
  
# This reactive computes the input for the cut off plot
  cutoff<-reactive({

      #call reactive trim
      # get dataframe conintaing the (average expression of a gene
      # in each treatment/condition group + maximum of this average value)
      matrix<-trim()[[3]]
      print("Inside module 2 line 269 cutof reactive")
      dat<-as.data.frame(matrix(NA, nrow = 200, ncol = 2))
      colnames(dat)<-c("x","y")

      #Loop to compute the number of genes that would remain in the dataset
      # for thresholds one to 200.
      for(i in 0:199)
      {
        dat[i+1,1]<-i #Threshold value (value for x-axis)
        dat[i+1,2]<-length(which(matrix$rowSum>=i)) #Get number of remaining genes after trimming (value for y-axis)
      }
      dat
  })
  
# Display cut off plot
  output$cut_off_plot <-renderPlotly({

      dat<-cutoff()
      # Plotly plot
      plot_ly(x=dat$x,y=dat$y,type = "scattergl",mode = 'markers',
              hoverinfo = 'text',
              text = ~paste('Cut-off:', dat$x,
                            '</br> Total gene counts:', dat$y)) %>%
        layout(
          xaxis = list(title = "Cut-off", gridcolor = "#bfbfbf", domain = c(0, 0.98)),
          yaxis = list(title = "Total gene counts", gridcolor = "#bfbfbf"))

  })
  
#Trim expression/count data using threshold if provided by user and update DESeq2 dataset 
dds.fc<-reactive({
    
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Preparing normalized table", value = 0)
    
    # Number of times we'll go through the loop to update progress bar
    n <- 3
    #Step 1: Retieve DESeq2 dataset object 
    dds<-trim()[[1]]
    
    #In order to visualize PCA we need to transform count data. Hence we perform rlogtransformation.
    #define variable to contain transformed data
    rld<-NULL
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    #Step 2: check if user has provided threshold to trim data, if not then skip trimming and proceed to normalization

      dds.fc<-dds
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      rld<-rlogTransformation(dds.fc)

    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    list(dds.fc,rld)
  })

#Once dataset has been trimmed by removing genes with low counts we proceed to the next step of pre-processing
#Perform normalization
#Get normalized data
normal <-reactive({

    dds <- dds.fc()[[1]]
    #Estimate size factors
    norm <- counts(dds, normalized=T)
    
    #Get normalized data

    if(!is.null(norm))
    {
      createAlert(session,"message","exampleAlert1", title="Message: To proceed to next step",
                  content = "Click on batch effect check button.", append=FALSE)
      
    }
    else closeAlert(session,"exampleAlert1")
    norm

})

# Display normalized table
output$normal <-
  DT::renderDataTable({
      DT::datatable(normal(),class = 'cell-border stripe',
                    selection = list(mode='single',target = 'row'),
                    extensions = list('Scroller'=NULL,'Buttons'=NULL),
                    options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 250,scroller = TRUE,dom = 'Bfrtip',
                                   buttons = list()))

  })

#display gene expression across conditon
observeEvent(input$normal_rows_selected,{
  #get which gene was clicked
  #get the row  number of the gene clicked
  selected_row <- input$normal_rows_selected
  print(selected_row)
  
  #get the normalized data
  full_data<-normal()

  #get the gene
  an_gene<-rownames(full_data)[selected_row]
  counts<-as.vector(full_data[selected_row,])
  if (as.numeric(filechoice())==3) condition <-cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice_module1())]) #already chosen for STAR
  else cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice)])

  library('data.table')
  df<-data.frame(counts,cond)
  colnames(df)<-c('count','condition')
  callModule(gene_count_module,"module",NULL,reactive({dds.fc()[[1]]}),reactive({an_gene}))
})

#download normalized data
output$downloadnorm <- downloadHandler(
  filename = function() {
    if(as.numeric(input$datachoice2==1)){paste("Normalized Data Table.xlsx")}
    else if (as.numeric(input$datachoice2==2)){paste("Normalized Data Table.csv")}
  },
  
  content = function(file) {
    if(as.numeric(input$datachoice2==1)){
      wb <- createWorkbook()
      addWorksheet(wb, sheetName = "Normalized table")
      writeData(wb = wb, sheet = 1, x = normal(), colNames = T, rowNames = T)
      saveWorkbook(wb, file)
      }
    else if (as.numeric(input$datachoice2==2)){write.csv(normal(), file)}
  }
)

conchoice_app<-reactive({
  if(as.numeric(filechoice() == 3)) conchoice_module1()
  else conchoice
})

  dds<-list(dds.fc=reactive({dds.fc()}),normal=reactive({normal()}),conchoice_app=reactive({conchoice_app()}))
  return(dds)
}
