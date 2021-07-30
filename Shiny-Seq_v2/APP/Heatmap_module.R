
heatmap_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    fluidRow(column(10, bsAlert("Anchor2"))),
    fluidRow(column(3,
                    selectInput(ns("Distance"), label = h5("Select distance method"), 
                                choices = list("Eucledian" = 1, "Manhattan" = 2),
                                selected = 1)
                    )),
    fluidRow(column(3,
                    selectInput(ns("Linkage"), label = h5("Select clustering method"), 
                                choices = list("Average" = 1,"complete"=2,"Ward.D2"=3,"Ward.D" = 4,"Singular"=5),
                                selected = 1)
                    )),
    uiOutput(ns("heat_comb")),
    downloadButton(ns('download_heatmap'), 'Download Plot'),
    fluidRow(column(8,imageOutput(ns("heatmap"))
                    )
             )
  )
}

heatmap_module<-function(input,output,session,heatmap_call,dds,batch,de_genes,combination,
                         CoCena, organism)
{
  print("line 29 heatmap")
  output$heat_comb <- renderUI({
    if(heatmap_call!="ANOVA" & heatmap_call != "TF")
    {
      result<-de_genes()
      combo<-combination()
      if(length(combo)>0)
      {
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Making plot", value = 0)
        
        # Number of times we'll go through the loop
        n <- 2
        
        rows<-length(combo)
        modules<-NULL
        res<-data.frame(matrix(NA, nrow = length(combo), ncol = 3))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
        if(!is.null(CoCena()))
        {
            cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
            entry<-c(combo, cluster_info$color)
            checklist<-list()
            for (i in seq_along(entry)) {
              checklist[[entry[[i]]]] = i
            }
            
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
            selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                        choices = checklist,selected = 1)
            
          
        }
        else{
          
          comb<-lapply(1:length(combo), function(i) {
            combo[[i]]
            
          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          
          # Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
          
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
          selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                      choices = checklist,selected = 1)
          
        }
        
      }
    } else if (heatmap_call == "TF"){
      if(organism() %in% c("Homo sapiens", "Mus musculus")){
        closeAlert(session,"Organism1")
      result<-de_genes()
      combo<-combination()
      if(length(combo)>0)
      {
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Making plot", value = 0)
        
        # Number of times we'll go through the loop
        n <- 2
        
        rows<-length(combo)
        modules<-NULL
        res<-data.frame(matrix(NA, nrow = length(combo), ncol = 3))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
        if(!is.null(CoCena()))
        {
            cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
            entry<-c(combo, cluster_info$color)
            checklist<-list()
            for (i in seq_along(entry)) {
              checklist[[entry[[i]]]] = i
            }
            
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
            selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                        choices = checklist,selected = 1)
            
          
        }
        else{
          
          comb<-lapply(1:length(combo), function(i) {
            combo[[i]]
            
          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          
          # Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
          
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
          selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                      choices = checklist,selected = 1)
          
        }
        
      }
      } else {
        createAlert(session, "Anchor2", "Organism1", title = "Error",
                    content = "This model is not available for the chosen organism")
      }
    }
    
  })

        heatmap<-reactive({
          
          dds.fc<-batch()$batch_data_for_DESeq()[[1]]
          rld<- assay(dds.fc)
          
          if(as.numeric(batch()$batch_choice())!=1)
          {
            rld<-batch()$batch_corrected()
          }
          
          if(heatmap_call=="ANOVA")
          {
            combo<-combination()
            num<-length(combo)
            heatmap_name<-"Heatmap of top 1000 most variable genes"
            heatmap_genes(heatmap_call,dds(),dds.fc,rld,NULL,
                          result,NULL,
                          num,heatmap_name,
                          input$Distance,input$Linkage)
          }
          else
          {
            result<-de_genes()
            combo<-combination()
            num<-length(combo)
            
            heatmap_name<-" "
            if(!is.null(CoCena()))
            {
                cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
                entry<-c(as.vector(combo), as.vector(cluster_info$color))
                
                heatmap_name<-paste("Heatmap of ",heatmap_call,entry[as.numeric(input$heat_choice)])
                h<-heatmap_genes(heatmap_call,NULL,dds.fc,rld,input$heat_choice,
                              result,
                              cluster_info,
                              num,heatmap_name,
                              input$Distance,input$Linkage)
                #insert validate
                shiny::validate(need(!is.null(h),'Heat map cannot be dislayed for one gene/trascription factor'))
                h
              
            }
            else
            {
              heatmap_name<-paste0("Heatmap of ",combo[[as.numeric(input$heat_choice)]])
              
              h<-heatmap_genes(heatmap_call,NULL,dds.fc,rld,input$heat_choice,
                            result,NULL,
                            num,heatmap_name,
                            input$Distance,input$Linkage)
              #insert validate
              shiny::validate(need(!is.null(h),'Heat map cannot be dislayed for one gene/trascription factor'))
              h
            }

          } 
          
        })

        #Display heatmap as .png with dimension 980 X 680
        output$heatmap <- renderImage({

          # A temp file to save the output. It will be deleted after renderImage
          # sends it, because deleteFile=TRUE.
          outfile <- tempfile(fileext='.png')
          # Generate a png
          if(!is.null(input$heat_choice)){
            # Create a Progress object
            progress <- shiny::Progress$new()
            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress$close())
            
            progress$set(message = "Making plot", value = 0)
            
            # Number of times we'll go through the loop
            n <- 2
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
            
            png(outfile, width=980, height=680)
            heatmap()
            dev.off()
          }
          else if(heatmap_call=='ANOVA')
          {
            png(outfile, width=980, height=680)
            heatmap()
            dev.off()
          }
          else
          {
            png(outfile, width=980, height=680)
            dev.off()
          }
          list(src = outfile,
               alt = "This is alternate text")

        },deleteFile=TRUE)
        
        
        #download heatmap
        output$download_heatmap <- downloadHandler(
          
          filename = function(){
            combo<-combination()
            choice<-as.numeric(input$heat_choice)
            heatmap_names<-""
            
            if(heatmap_call=="ANOVA")  heatmap_name<-"Heatmap of top 1000 most variable genes"
            else heatmap_name<-paste("Heatmap of ",heatmap_call,combo[[choice]])
            paste(heatmap_name,'.pdf')
            
          },
          content = function(file) {
            
            png(file,width = 980, height = 680)
             heatmap()
            dev.off()
            
          })   
    
return(list(
  input_Distance=reactive({input$Distance}),
  input_Linkage=reactive({input$Linkage})
))         
        
}
  