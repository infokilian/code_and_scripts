
source("gene_count_module.R")
ANOVA_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    br(),
    radioButtons(ns("lfc"), "Show log2 foldchange?:", choices = c("Yes", "No"),selected = "No"),
    br(),
    actionButton(ns("ok4"), "Start"),
    br(),
    br(),
    fluidRow(
      column(1,
             selectInput(ns("datachoice6")  ,label = h5("Select Data Type"), 
                         choices = list("Excel" = 1, "CSV" = 2),
                         selected = 1)),
      column(1, 
             br(),
             br(),
             
             downloadButton(ns('download_ANOVA_Table'), 'Download ANOVA Table'))),
    fluidRow(column(5,bsAlert("DE_missing"))),
    DT::dataTableOutput(ns("ANOVA")),
    gene_count_module_UI(ns("module5"))
  )
}
ANOVA_module<-function(input,output,session,
                       batch_choice,dds,design,combination,
                       conchoice,DE_genes,normal,batch_corrected,ok3_de)
{
  
  observeEvent(ok3_de(),{
  if(ok3_de()>0) {
    closeAlert(session, "message1")
  } else {
    createAlert(session,"DE_missing",alertId = "message1",
                title = "Error",
                content = "Please run the Differential Expression Analysis before computing the ANOVA")
  }
})

    print("inside anova module")
    
    #Perform ANOVA 
    anova<-reactive({
      res<-NULL
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Preparing ANOVA table", value = 0)
      
      # Number of times we'll go through the loop
      n <- 2
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      if(as.numeric(batch_choice())==2)
      {
        
        res<-DESeq(dds(),test = "LRT",reduced = design())
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      else if(as.numeric(batch_choice())==3)
      {
        res<-DESeq(dds(),test = "LRT",reduced = design())
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      else 
      {
        res<-DESeq(dds(),test = "LRT",reduced = design())#,parallel = TRUE)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      res
    })
    
    #prepare ANOVA table for display(This reactive merges the anova table with DE analysis results)
    anova_table<-reactive({
      
      dds_anova<-anova()
      res<-results(dds_anova)

      #Get base mean of each condition
      base_mean<-sapply( levels(dds_anova$condition), function(lvl) rowMeans( counts(dds_anova,normalized=TRUE)[,dds_anova$condition == lvl] ) )

      #Merge the basemeans with the anova table
      anova<-merge(res,base_mean,by=0,all=TRUE)
      for(i in 8:length(colnames(anova)))
      {
        colnames(anova)[i]<-paste0(colnames(anova)[i],'_mean')
      }

      genes<-anova[,1]
      anova<-anova[-1]
      rownames(anova)<-genes
      colnames(anova)[1]<-"Overall mean"
      result<-DE_genes()
      combo<-combination()
      num <-length(combo)
      l<-list()
      l[[1]]<-anova[,-3]
      for(i in 1:num)
      {
        res<-as.data.frame(result()[[i]][4])
        ######change####
        a_tab<-res[,-4]
        c<-colnames(a_tab)
        res <- a_tab[,c(c[2],c[3],c[5],c[6],c[1],c[4])]
        ####################
        
        colnames(res)[5]<-"Overall mean"
        for (j in 1:length(colnames(res)))
        {
          colnames(res)[j]<-paste0(combo[[i]],' ',colnames(res)[j])
        }
        
        l[[length(l)+1]]<-res
        
      }
      
      #Condense results of ANOVA and DE analysis into one data
      anova<-transform(Reduce(merge, lapply(l, function(x) data.frame(rn = rownames(x),x))), row.names=rn, rn=NULL)
    })
    
    
 observeEvent(input$ok4,{    
    #Display ANOVA table
    output$ANOVA <- DT::renderDataTable({
      a_tab<-anova_table()[,-c(2,3)]
      ########change#######
      cond<-unique(colData(dds())[,as.numeric(conchoice())])
      c<-colnames(a_tab)
      temp<-as.vector(c[4:(3+length(cond))])
      temp2<-as.vector(c[(length(cond)+4):length(c)])
      ######################
      a_tab=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
      anova<-NULL
      if(input$lfc =="No")
      {
        names<-colnames(a_tab)
        idx<-grep("log2FoldChange",names)
        anova=a_tab[,-idx]
      }
      else anova=a_tab
      
      DT::datatable(anova,class = 'cell-border stripe',
                    selection = list(mode='single',target = 'row'),
                    extensions = list('Scroller'=NULL,'Buttons'=NULL),
                    options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 300,scroller = TRUE,dom = 'Bfrtip',
                                   buttons = list()
                    )
      )
      
    })
    #display boxplot when row(gene) in ANOVA table is clicked
    observeEvent(input$ANOVA_rows_selected,{
      selected <- input$ANOVA_rows_selected
      full_data<-normal()

      #get the gene
      temp<-strsplit(rownames(anova_table())," ")
      an_gene<-temp[[selected]]
      library('data.table')
      counts<-as.vector(full_data[which(an_gene %in% rownames(full_data)),])
      cond<-as.vector(colData(dds())[,as.numeric(conchoice())])
      df<-data.frame(counts,cond)
      colnames(df)<-c('count','condition')
      batch<-batch_choice()
      if(as.numeric(batch)==1) callModule(gene_count_module,"module5",NULL,reactive({dds()}),reactive({an_gene}))
      else
        {
          idx<-which(rownames(batch_corrected()) %in% an_gene)
          print(idx)
          callModule(gene_count_module,"module5",batch_corrected()[idx,],reactive({dds()}),reactive({an_gene}))
        }
    })
    
    #Download ANOVA table
    output$download_ANOVA_Table <- downloadHandler(
      
      filename = function() { 
        if(as.numeric(input$datachoice6==1)){paste('ANOVA_table.xlsx')}
        else if(as.numeric(input$datachoice6==2)){paste('ANOVA_table.csv')}
        },
      content = function(file) {
        
      if(as.numeric(input$datachoice6==1)){
        anova<-anova_table()
        a_tab<-anova_table()[,-c(2,3)]
        ########reordering columns in anova table for display#######
        cond<-unique(colData(dds())[,as.numeric(conchoice())])
        c<-colnames(a_tab)
        temp<-as.vector(c[4:(3+length(cond))])
        temp2<-as.vector(c[(length(cond)+4):length(c)])
        ######################
        anova=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
        M <- as.matrix(anova)
        wb <- createWorkbook()
        addWorksheet(wb, sheetName = "ANOVA")
        writeData(wb = wb, sheet = 1, x = M, colNames = T, rowNames = T)
        saveWorkbook(wb, file)
      }
      else {
        anova<-anova_table()
        a_tab<-anova_table()[,-c(2,3)]
        ########reordering columns in anova table for display#######
        cond<-unique(colData(dds())[,as.numeric(conchoice())])
        c<-colnames(a_tab)
        temp<-as.vector(c[4:(3+length(cond))])
        temp2<-as.vector(c[(length(cond)+4):length(c)])
        ######################
        anova=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
        write.csv(anova, file)
      }
        }
    )

  }) 
  return(list(
    anova_table=reactive({anova_table()}),
    dds=reactive({anova()})
  )
  ) 
  
}
