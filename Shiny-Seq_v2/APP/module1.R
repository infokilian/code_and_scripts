Module_Raw_data_UI <- function(id)

{
  ns <- NS(id)
  tagList(
    useShinyjs(),
    br(),
    # Select organism
    fluidRow(column(2,
                    br(),
                    br(),
                    selectInput(ns("organism"), label = h5("Select organism"), 
                                choices = list("", "Homo sapiens", 
                                               "Mus musculus", 
                                               "Anopheles gambiae",
                                               "Arabidopsis thaliana",
                                               "Bos taurus",
                                               "Caenorhabditis elegans",
                                               "Canis familiaris",
                                               "Danio rerio",
                                               "Drosophila melanogaster",
                                               "Escherichia coli K12" = "K-12",
                                               "Escherichia coli Sakai" = "Sakai",
                                               "Gallus gallus",
                                               "Macaca mulatta",
                                               "Pan troglodytes",
                                               "Plasmodium falciparum",
                                               "Rattus norvegicus",
                                               "Sus scrofa",
                                               "Xenopus laevis",
                                               "Others"
                                               ),
                                selected = NULL))),

    # Select reference genome
    fluidRow(column(2, 
                    selectInput(ns("refgenome"), label = h5("Select reference genome"),
                                choices = list("", "BioMart" = 1, "Own upload" = 2),
                                selected = NULL)),                    
             column(1,
                    br(),
                    br(),
                    actionButton(ns("help1"), label = "Help", icon = icon("question-circle")))),

    # Choose starting point
    uiOutput(ns("inp_ch")),
    
    # Upload reference genome
    fluidRow(column(4, uiOutput(ns("file_ref_genome")))),
    
    #Upload STAR-aligned data
    fluidRow(column(12, uiOutput(ns("file4")))),

    # Upload Kallisto data
    fluidRow(column(4, uiOutput(ns("file3")))),

    # Upload expression table
    fluidRow(column(4, uiOutput(ns("file1")))),

    # Upload annotation table
    fluidRow(column(12, uiOutput(ns("file2")))),
    
    # if STAR, specify condition and design in Raw data
    fluidRow(column(4, uiOutput(ns("condition_module1")))),
    fluidRow(column(4, uiOutput(ns("design_module1")))),
    

    bsModal(ns("modalhelp"), "Selection of a reference genome", ns("help1"),  
            helpText('When using Kallisto data, a gene annotation file is needed, which
                     links transcript identifiers to gene names.\n
                     Either select the BioMart database to construct a gene annotation file 
                     or upload your own file.\n
                     In case of uploading your own file, it must be a .csv or .tsv file 
                     containing a column called "TXNAME" which holds the ensembl transcript id versions 
                     and a column called "SYMBOL" containing the external gene names / gene symbols.'))
  

  )
}

Module_Raw_data_Input <- function(input, output, session)
{

  # Aim: Display count/expression table and annotation table
    # Option 1: Generate count/expression table from Kallisto zip file and upload annotation table
    # Option 2: Upload count/expression table and annotation table
    # Option 3: Generate count/expression table from STAR zip file and upload annotation table

  print("module 1 (Raw Data) line 62 - start of server - ok")

  output$inp_ch <- renderUI({
    if(input$organism != "")
      fluidRow(
        column(2,
               radioButtons(session$ns("filechoice"), "Choose starting point",
                            choices = list("Kallisto" = 1, "Count table" = 2, "STAR" = 3), selected = 1)),
        column(2,
               radioButtons(session$ns("sep"), "Separator",
                            choices = c(Comma = ",",
                                        Space = " ",
                                        Tab = "\t")))
      )
  })

  input_file <- function(in_file)
  {
    if(is.null(in_file))
    {return(NULL)}
    else
      return(read.table(in_file$datapath, header = TRUE, sep = input$sep, check.names = FALSE, quote = "\""))
  }

  # Option 1: Generate count/expression table from Kallisto zip file and upload annotation table

    # Upload Kallisto zip file to compute count/expression table

  output$file3 <- renderUI({
    req(input$filechoice)
    if(as.numeric(input$filechoice) == 1 && input$organism != "")
      {fileInput(session$ns("kallistozipfile"), label = h5("Upload Kallisto zip file"),
                 accept = ".zip")}
  })

  observe({
    req(input$kallistozipfile)
    if(is.null(input$kallistozipfile))
      return(NULL)
    else
      return(unzip(input$kallistozipfile$datapath,
                   exdir = paste0("./unzipped_kallistofiles/", session$token)))
  })

  output$file_ref_genome <- renderUI({
    if(as.numeric(input$refgenome == 2)){
      fileInput(session$ns("refgenome_path"), label = h5('Choose reference genome file to upload'),
                accept = c('.csv',
                           '.tsv',
                           'text/comma-separated-values',
                           'text/csv',
                           'text/plain',
                           'text/tab-separated-values')
                )
      }
    })
  
  orgDB <- reactive({
    
    db <- NULL
    if(input$organism == "Homo sapiens") db <- org.Hs.eg.db
    else if(input$organism == "Mus musculus") db <- org.Mm.eg.db
    else if(input$organism == "Anopheles gambiae") db <- org.Ag.eg.db
    else if(input$organism == "Arabidopsis thaliana") db <- org.At.tair.db
    else if(input$organism == "Bos taurus") db <- org.Bt.eg.db
    else if(input$organism == "Caenorhabditis elegans") db <- org.Ce.eg.db
    else if(input$organism == "Canis familiaris") db <- org.Cf.eg.db
    else if(input$organism == "Danio rerio") db <- org.Dr.eg.db
    else if(input$organism == "Drosophila melanogaster") db <- org.Dm.eg.db
    else if(input$organism == "K-12") db <- org.EcK12.eg.db
    else if(input$organism == "Sakai") db <- org.EcSakai.eg.db
    else if(input$organism == "Gallus gallus") db <- org.Gg.eg.db
    else if(input$organism == "Macaca mulatta") db <- org.Mmu.eg.db
    else if(input$organism == "Pan troglodytes") db <- org.Pt.eg.db
    else if(input$organism == "Plasmodium falciparum") db <- org.Pf.plasmo.db
    else if(input$organism == "Rattus norvegicus") db <- org.Rt.eg.db
    else if(input$organism == "Sus scrofa") db <- org.Ss.eg.db
    else if(input$organism == "Xenopus laevis") db <- org.Xl.eg.db #02/01/2021 laevis not available
    db
  })
   
  dataset <- reactive({

    org <- NULL
    if(input$organism == "Homo sapiens") org <- "hsapiens_gene_ensembl"
    else if(input$organism == "Mus musculus") org <- "mmusculus_gene_ensembl"
    else if(input$organism == "Bos taurus") org <- "btaurus_gene_ensembl"
    else if(input$organism == "Caenorhabditis elegans") org <- "celegans_gene_ensembl"
    else if(input$organism == "Danio rerio") org <- "drerio_gene_ensembl"
    else if(input$organism == "Drosophila melanogaster") org <- "dmelanogaster_gene_ensembl"
    else if(input$organism == "Gallus gallus") org <- "ggallus_gene_ensembl"
    else if(input$organism == "Macaca mulatta") org <- "mmulatta_gene_ensembl"
    else if(input$organism == "Pan troglodytes") org <- "ptroglodytes_gene_ensembl"
    else if(input$organism == "Rattus norvegicus") org <- "rnorvegicus_gene_ensembl"
    else if(input$organism == "Sus scrofa") org <- "sscrofa_gene_ensembl"
    else if(input$organism == "Xenopus laevis") org <- "xtropicalis_gene_ensembl" #02/01/2021 laevis not available
    org
  })
  
  # This reactive prepares input for the reactive txi which generates count table from Kallisto files
  t2g <- reactive({
    if(!exists("mart") & !is.null(dataset()))
      {
      mart <- biomaRt::useMart(host = "https://nov2020.archive.ensembl.org",
                               biomart = "ENSEMBL_MART_ENSEMBL",
                               dataset = dataset())
      t2g <- biomaRt::getBM(
        attributes = c("ensembl_transcript_id_version",
                       "ensembl_gene_id",
                       "ensembl_gene_id_version",
                       "external_gene_name"),
        mart = mart)
      t2g
      }
    })

    # Create the count matrix out of the Kallisto abundance table
  txi <-
    reactive({
      if(as.numeric(input$filechoice) == 1 && input$kallistozipfile != "")
        {
        # Create a progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Processing Kallisto files", value = 0)

        print("Directory of the uploaded Kallisto zip file:")
        print(dir(paste0("./unzipped_kallistofiles/", session$token)))

        # Increment the progress bar, and update the detail text
        progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
        # Pause for 0.1 seconds to simulate a long computation
        Sys.sleep(0.1)

        if(as.numeric(input$refgenome == 1))
          {
          print("tx2gene")
          print(t2g()[, c("ensembl_transcript_id_version",
                          "external_gene_name")])
          d <- tximport(sapply(dir(paste0("./unzipped_kallistofiles/", session$token)),
                               function(id)
                                 {
                                 print(list.files(file.path("./unzipped_kallistofiles", session$token, id),
                                                  "*.h5", full.names = TRUE))
                                 }
                               ),
                        "kallisto",
                        tx2gene = t2g()[, c("ensembl_transcript_id_version",
                                            "external_gene_name")]
                        )
          }
        else
          {
          t2x <- read.delim(file = input$refgenome_path$datapath, sep = ",")
          t2x <- t2x[, c("TXNAME", "SYMBOL")]
          d <- tximport(sapply(dir(paste0("./unzipped_kallistofiles/", session$token)),
                               function(id)
                                 {
                                 print(list.files(file.path("./unzipped_kallistofiles", session$token, id),
                                                  "*.h5", full.names = TRUE))
                                 }
                               ),
                        "kallisto",
                        tx2gene = t2x
                        )
          }

        # Increment the progress bar, and update the detail text
        progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
        # Pause for 0.1 seconds to simulate a long computation
        Sys.sleep(0.1)

        print("module 1 (Raw Data) line 189 - Count matrix sucessfully created out of the Kallisto abundance table")
        d
        }
      })

  # Option 2: Upload count/expression table and annotation table

    # Input file containing count/expression table from user
  output$file1 <- renderUI({
    req(input$filechoice)
    if(as.numeric(input$filechoice) == 2){
      fileInput(session$ns("file1"), label = h5('Choose Expression file to upload'),
                accept = c('.csv',
                           '.tsv',
                           'text/comma-separated-values',
                           'text/csv',
                           'text/plain',
                           'text/tab-separated-values')
                )
      }
    })

  # Option 3: Generate count/expression table from STAR zip file and upload annotation table 
  # Upload STAR zip file to compute count/expression table
  output$file4 <- renderUI({
    req(input$filechoice)
    if(as.numeric(input$filechoice) == 3 && input$organism != "")
    {
      fluidRow(column(2,fileInput(session$ns("STARzipfile"), label = h5("Upload STAR zip file"),
                        accept = ".zip")),
               column(1,
                      br(),
                      br(),
                      actionButton(session$ns("help2"), "Help", icon = icon("question-circle"))),
               bsModal(session$ns("modalhelp2"), "Structure of the annotation file", session$ns("help2"),
                       helpText(HTML(paste("In order to ensure the correct read in of the STAR aligned data, the uploaded .zip-file should include only
                                the reads-per-gene tables","\n", h4("Do not store the reads-per-gene tables in a sub-folder!")))))
               )
               
    }
  })
  
  #if STAR is chosen choose design in Raw data
  # creates a checkbox widget to select the important condition
  output$condition_module1 <-
    renderUI({
      req(input$file2)
      req(input$filechoice)
      if(input$filechoice == 3){
        checklist = list()
        for (i in seq_along(colnames(input_file(input$file2)))[-1]) {
          checklist[[colnames(input_file(input$file2))[[i]]]] = i
        }
        radioButtons(session$ns("conchoice_module1"), "Choose the conditon", checklist,selected = 2)
      }
    })
  
  # creates a checkbox widget to select the variables to be included in design
  output$design_module1 <-
    renderUI({
      req(input$file2)
      req(input$filechoice)
      req(input$conchoice)
      if(input$filechoice==3){
        checklist = list()
        for (i in seq_along(colnames(input_file(input$file2)))[-1]) {
          # disply all variables except condition variable
          if(i!=as.numeric(input$conchoice)) checklist[[colnames(input_file(input$file2))[[i]]]] = i
        }
        if(length(checklist)!=0){
          checkboxGroupInput(session$ns("designchoice_module1"), "Choose the variables to be included in design", checklist)
        }
      }
    })
  
  observe({
    req(input$STARzipfile)
    if(is.null(input$STARzipfile))
      return(NULL)
    else
      return(unzip(input$STARzipfile$datapath,
                   exdir = paste0("./unzipped_STARfiles/", session$token)))
  })
  
  dds_txi <- reactive({

    if(as.numeric(input$filechoice) == 3 && input$STARzipfile != ""){
    req(input$file2)
      
      # Create a progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Processing STAR files", value = 0)
      
    directory <- file.path("./unzipped_STARfiles/", session$token) #set directory of the STAR count data
    sample_table <- input_file(input$file2)
    
    condition<-as.numeric(input$conchoice_module1)
    design<-as.numeric(input$designchoice_module1)
    colnames(sample_table)[condition]<-"condition"
    #update design formula
    d<-''
    for (i in design)
    {
      d<-paste(d,colnames(sample_table)[i]," + ",sep=" ")
    }
    d<-paste('~',' ',d,colnames(sample_table)[condition],sep = " ")
    
    sample_table$Files <- paste(sample_table[,1], ".ReadsPerGene.out.tab", sep = "") #create a vector consisting of the respective file names
    sample_table <- sample_table[,c(1,ncol(sample_table),3:ncol(sample_table)-1)]
    dds_txi <- as.data.frame(counts(DESeqDataSetFromHTSeqCount(sampleTable = sample_table, #create the count table
                                          directory = directory,
                                          design = formula(d))))
    dds_txi <- dds_txi[5:nrow(dds_txi),] #delete unnecessary rows
    
    # Increment the progress bar, and update the detail text
    progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
    # Pause for 0.1 seconds to simulate a long computation
    Sys.sleep(0.1)
    
    if(as.numeric(input$refgenome == 1))
    {
      #annotate the count table with the respective gene symbols
      dds_txi$ensembl_id <- rownames(dds_txi)
      dds_txi_anno <- merge(dds_txi,
                            t2g()[, c("ensembl_gene_id_version",
                                      "external_gene_name")],
                         by.x = "ensembl_id",
                         by.y = "ensembl_gene_id_version")
      print("dds_txi")
      print(dds_txi_anno)
      dds_txi_anno <- dds_txi_anno[!duplicated(dds_txi_anno$external_gene_name),]
      rownames(dds_txi_anno) <- dds_txi_anno$external_gene_name
      dds_txi_anno$ensembl_id <- NULL
      dds_txi_anno$external_gene_name <- NULL
      dds_txi_anno <- dds_txi_anno[order(rownames(dds_txi_anno), decreasing = F),]
    } else {
      t2x <- read.delim(file = input$refgenome_path$datapath, sep = ",")
      t2x <- t2x[, c("TXNAME", "SYMBOL")]
      dds_txi$txname <- rownames(dds_txi)
      dds_txi_anno <- merge(dds_txi,
                            t2x,
                            by.x = "txname",
                            by.y = "TXNAME")
      dds_txi_anno <- dds_txi_anno[!duplicated(dds_txi_anno$SYMBOL),]
      rownames(dds_txi_anno) <- dds_txi_anno$SYMBOL
      dds_txi_anno$txname <- NULL
      dds_txi_anno$TXNAME <- NULL
      dds_txi_anno <- dds_txi_anno[order(rownames(dds_txi_anno), decreasing = F),]
    }
    # Increment the progress bar, and update the detail text
    progress$inc(1/2, detail = paste("Doing part", 2,"/",2))
    # Pause for 0.1 seconds to simulate a long computation
    Sys.sleep(0.1)
    
    dds_txi_anno
    }
  })

  #Input file containing annotation file as input from the user
                 output$file2<-renderUI({
 
                   req(input$filechoice)
                   if(as.numeric(input$filechoice)==2)
                   {
                     req(input$file1)
                     fluidRow(column(2, fileInput(session$ns("file2"), label = h5('Choose Annotation file to upload'),
                                                  accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv'
                                                  )
                                    )),
                              column(1,
                                     br(),
                                     br(),
                                     actionButton(session$ns("help2"), "Help", icon = icon("question-circle"))),
                              bsModal(session$ns("modalhelp2"), "Structure of the annotation file", session$ns("help2"),  
                                      helpText(HTML(paste("The annotation file contains additional metainformation on the samples of the uploaded expression file.
                                               In order to ensure an accurate assignment of the metainformation to the correct sample, the annotation file
                                               should be structured as follows:",
                                               "1. row: The first row contains the name of the variables (starting with: SampleID).",
                                               "1. column: SampleID - This column contains the names of the samples which were also used in the expression file.",
                                               "2.- i. column: NameOfVariable - Theses columns contain the metainformation for each sample regarding the variables defined in the first row.", 
                                               sep="<br/>"))),
                                      fluidRow(column(5,img(src="anno_table.PNG"))))
                              )
                   }
                   else if( as.numeric(input$filechoice)==1)
                   {
                     req(input$kallistozipfile)
                     fluidRow(column(2, fileInput(session$ns("file2"), label = h5('Choose Annotation file to upload'),
                                                  accept = c(
                                                    'text/csv',
                                                    'text/comma-separated-values',
                                                    'text/tab-separated-values',
                                                    'text/plain',
                                                    '.csv',
                                                    '.tsv'
                                                  )
                     )),
                     column(1,
                            br(),
                            br(),
                            actionButton(session$ns("help2"), "Help", icon = icon("question-circle"))),
                     bsModal(session$ns("modalhelp2"), "Structure of the annotation file", session$ns("help2"),  
                             helpText(HTML(paste("The annotation file contains additional metainformation on the samples of the uploaded expression file.
                                               In order to ensure an accurate assignment of the metainformation to the correct sample, the annotation file
                                               should be structured as follows:",
                                               "1. row: The first row contains the name of the variables (starting with: SampleID).",
                                               "1. column: SampleID - This column contains the names of the samples which were also used in the expression file.",
                                               "2.- i. column: NameOfVariable - Theses columns contain the metainformation for each sample regarding the variables defined in the first row.", 
                                               sep="<br/>"))),
                             fluidRow(column(5,img(src="anno_table.PNG"))))
                     )
                   }
                   else if( as.numeric(input$filechoice)==3)
                   {
                     req(input$STARzipfile)
                     fluidRow(column(2, fileInput(session$ns("file2"), label = h5('Choose Annotation file to upload'),
                                                  accept = c(
                                                    'text/csv',
                                                    'text/comma-separated-values',
                                                    'text/tab-separated-values',
                                                    'text/plain',
                                                    '.csv',
                                                    '.tsv'
                                                  )
                     )),
                     column(1,
                            br(),
                            br(),
                            actionButton(session$ns("help2"), "Help", icon = icon("question-circle"))),
                     bsModal(session$ns("modalhelp2"), "Structure of the annotation file", session$ns("help2"),  
                             helpText(HTML(paste("The annotation file contains additional metainformation on the samples of the uploaded expression file.
                                               In order to ensure an accurate assignment of the metainformation to the correct sample, the annotation file
                                               should be structured as follows:",
                                                 "1. row: The first row contains the name of the variables (starting with: SampleID).",
                                                 "1. column: SampleID - This column contains the names of the samples which were also used in the expression file.",
                                                 "2.- i. column: NameOfVariable - Theses columns contain the metainformation for each sample regarding the variables defined in the first row.", 
                                                 sep="<br/>"))),
                             fluidRow(column(5,img(src="anno_table.PNG"))))
                     )
                   }
                 })

  #compose count /expression table obtained as input to be displayed
  #edata reactive contains the count/expression table
  edata<-reactive({
    #setting up expression data
    #Get expression table from file
    
    if(as.numeric(input$filechoice)==2)
    {
      input_exp<-input_file(input$file1)
      if(!is.null(input_exp))
      {
        data<-input_exp[,-1]
        col_names<-colnames(data)
        #convert entries to integers
        data=as.matrix(data)
        storage.mode(data)="integer"
        data= data.frame(data)
        rownames(data)=as.character(input_exp[,1])
        colnames(data)<-col_names
        data<-data
      }
    }
    else if(as.numeric(input$filechoice)==1)
    {
      dat<-as.data.frame(txi()$counts)
      sample_names<-colnames(dat)
      #convert entries to integers
      data=as.matrix(dat)
      storage.mode(data)="integer"
      data= data.frame(data)
      rownames(data)=rownames(dat)
      colnames(data)<-sample_names

      data <- data

    }else if(as.numeric(input$filechoice)==3)
    {
      dat<-as.data.frame(dds_txi())
      sample_names<-colnames(dat)
      #convert entries to integers
      data=as.matrix(dat)
      storage.mode(data)="integer"
      data= data.frame(data)
      rownames(data)=rownames(dat)
      colnames(data)<-sample_names
      
      data <- data
      
    }

  })
  
  edata2 <- reactive({
    if(as.numeric(input$filechoice)==2)
    {
      input_exp<-input_file(input$file1)
      if(!is.null(input_exp))
      {
        data<-input_exp[,-1]
        col_names<-colnames(data)
        #convert entries to integers
        data=as.matrix(data)
        storage.mode(data)="integer"
        data= data.frame(data)
        rownames(data)=as.character(input_exp[,1])
        colnames(data)<-col_names
        data<-data
      }
    }
    else if(as.numeric(input$filechoice)==1)
    {
      data <- txi()
    } 
    else if(as.numeric(input$filechoice)==3)
    {
      data <- dds_txi()
    }
  })

  #Compose annotation data from file
  #pdata contains annotation file
  pData<-reactive({
    #setting up annotation data
    #Get annotation table from file
    pheno<-input_file(input$file2)

    #Get expression table
    data<-edata()

    #Expression table and annotation table should not be null
    if(!is.null(data) && !is.null(pheno))
    {

      #Assuming that the first column of annotation table is sample id ,extract all sample IDs
      sample_id = pheno[,1]
      #Get all sample IDs from expression table(sample ID refer to column names of expression table)
      exp_sample_id = colnames(data)
      #Check if all sample ID in expression table are present in the annotation table
      if (all(exp_sample_id %in% sample_id))
      {
        #set all variables of annotation table as factors
        col<-1:ncol(pheno)
        for (i in col)
        {
          pheno[,i]<-as.factor(pheno[,i])
        }
        #Remove those sample IDs that are present in the expression table but absent in the annotation table from the annotation table
        idx<-NULL
        if(!all(sample_id %in% exp_sample_id))
        {
          idx<-which(!(sample_id %in% exp_sample_id))
          pheno <- pheno[-idx, ]
        }

        pheno<-pheno

      }
      else 
        {
          pheno<-NULL
        }
      pheno
    }
  })
out<-list(file2=reactive({input$file2}),
          edata=reactive({edata()}),
          pData=reactive({pData()}),
          organism=reactive({input$organism}),
          edata2=reactive({edata2()}),
          refgenome=reactive({input$refgenome}),
          dataset=reactive({dataset()}),
          orgDB=reactive({orgDB()}),
          filechoice=reactive({input$filechoice}),
          condition_module1=reactive({input$conchoice_module1}),
          design_module1=reactive({input$designchoice_module1}))
return(out)
}

onStop(fun = function() {
  if(input$filechoice==1){
  files_info <- file.info(list.files(path = "./unzipped_kallistofiles",
                                    full.names = TRUE),
                         extra_cols = FALSE)
  delete_files <- files_info[difftime(Sys.time(), 
                                     files_info[,"mtime"],
                                     units = "days") > 1 , 1:4]
  unlink((row.names(delete_files)), recursive = TRUE)
  } else if (input$filechoice==3){
  files_info <- file.info(list.files(path = "./unzipped_STARfiles",
                                     full.names = TRUE),
                          extra_cols = FALSE)
  delete_files <- files_info[difftime(Sys.time(), 
                                      files_info[,"mtime"],
                                      units = "days") > 1 , 1:4]
  unlink((row.names(delete_files)), recursive = TRUE)
  }
})