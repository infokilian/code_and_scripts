source("functions.R")
#Module PCA
#UI component
#Server component: Calls Pcaplot function from fucntions.R

PCA_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    #User input: PCA type: 2D/3D
    radioButtons(ns("plotType"), "PCA Type:", choices = c("2D", "3D")),
    #User input: Annotation variables to color samples in PCA
    uiOutput(ns("annotation_pca")),
    #User input: Compute PCA for top 500 most variable genes/All genes
    selectInput(ns("top"), label = h5("Should PCA be computed for:"), 
                choices = list("Top 500 genes" = 1, "All" = 2),
                selected = 2),
    uiOutput(ns("point_size")),
    #Download button: PCA gets downloaded on user click
    fluidRow(column(3," "),column(6,downloadButton(ns('download_pca_svg'), 'Download plot as svg'))),
    #Output: Display PCA plot
    fluidRow(column(1," "),
             column(6,plotlyOutput(ns("pca"))))
  )
}

PCA_module<-function(input,output,session,data,resize_factor=NULL)
{
  output$point_size<-renderUI({
    req(input$plotType)
    if(input$plotType=="2D")
    {
      textInput(session$ns("pca_2d_plot_point_size"),label="Set point size in plot",value = "5")
    }
  })
     #Creates a checkbox widget containing annotation variables to color samples in PCA
     output$annotation_pca<-renderUI({
                            checklist = list()
                            for (i in seq_along(colnames(data[[2]]))) { 
                                  checklist[[colnames(data[[2]])[[i]]]] = i
                                 }
                             radioButtons(session$ns("annotation_pca"), "Choose the annotation", checklist)
                           })
  
     #Generate PCA plot
      output$pca<-renderPlotly({
                  point_size=1
                  req(input$pca_2d_plot_point_size)
                  if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
                  pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
                          plotType = input$plotType, resize_factor = resize_factor,
                          point_size = point_size)[[1]]
  
                  })
      
    #Download pca in svg format
      output$download_pca_svg <- downloadHandler(
  
                                 filename = function(){
                                   
                                   if (input$plotType == "2D")
                                   {
                                     if(input$top==2)  paste('2D PCA of Raw data all genes.svg')
                                     else  paste('2D PCA of Raw data top 500 genes.svg')
                                   }
                                   else if(input$plotType == "3D")
                                   {
                                     if(input$top==2)  paste('3D PCA of Raw data all genes.pdf')
                                     else  paste('3D PCA of Raw data top 500 genes.pdf')
                                   }
                                 },
                                 
                                 content = function(file) {
                                   point_size=1
                                   req(input$pca_2d_plot_point_size)
                                   if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
                                   
                                   p<-pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
                                              plotType = input$plotType, point_size = point_size)[[2]]
                                   ggsave(file,p, width = 22, height = 15, units = "cm")
                                 }
                                 )
      return(reactive({input$top}))
}

PCA_module2<-function(input,output,session,data,resize_factor=NULL)
{
  output$point_size<-renderUI({
    req(input$plotType)
    if(input$plotType=="2D")
    {
      textInput(session$ns("pca_2d_plot_point_size"),label="Set point size in plot",value = "5")
    }
  })
  #Creates a checkbox widget containing annotation variables to color samples in PCA
  output$annotation_pca<-renderUI({
    checklist = list()
    for (i in seq_along(colnames(data[[2]]))) { 
      checklist[[colnames(data[[2]])[[i]]]] = i
    }
    radioButtons(session$ns("annotation_pca"), "Choose the annotation", checklist)
  })
  
  #Generate PCA plot
  output$pca<-renderPlotly({
    point_size=1
    req(input$pca_2d_plot_point_size)
    if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
    pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
            plotType = input$plotType, resize_factor = resize_factor,
            point_size = point_size)[[1]]
    
  })
  #Download pca in svg format
  output$download_pca_svg <- downloadHandler(
    
    filename = function(){
      
      if (input$plotType == "2D")
      {
        if(input$top==2)  paste('2D PCA of Normalized data all genes.svg')
        else  paste('2D PCA of Normalized data top 500 genes.svg')
      }
      else if(input$plotType == "3D")
      {
        if(input$top==2)  paste('3D PCA of Normalized data all genes.pdf')
        else  paste('3D PCA of Normalized data top 500 genes.pdf')
      }
    },
    
    content = function(file) {
      point_size=1
      req(input$pca_2d_plot_point_size)
      if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
      
      p<-pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
                 plotType = input$plotType, point_size = point_size)[[2]]
      ggsave(file,p, width = 22, height = 15, units = "cm")
    }
  )
  return(reactive({input$top}))
}

PCA_module3<-function(input,output,session,data,resize_factor=NULL)
{
  output$point_size<-renderUI({
    req(input$plotType)
    if(input$plotType=="2D")
    {
      textInput(session$ns("pca_2d_plot_point_size"),label="Set point size in plot",value = "5")
    }
  })
  #Creates a checkbox widget containing annotation variables to color samples in PCA
  output$annotation_pca<-renderUI({
    checklist = list()
    for (i in seq_along(colnames(data[[2]]))) { 
      checklist[[colnames(data[[2]])[[i]]]] = i
    }
    radioButtons(session$ns("annotation_pca"), "Choose the annotation", checklist)

  })
  
  #Generate PCA plot
  output$pca<-renderPlotly({
    point_size=1
    req(input$pca_2d_plot_point_size)
    if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
    pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
            plotType = input$plotType, resize_factor = resize_factor,
            point_size = point_size)[[1]]
    
  })
  
  #Download pca in svg format
  output$download_pca_svg <- downloadHandler(
    
    filename = function(){
      
      if (input$plotType == "2D")
      {
        if(input$top==2)  paste('2D PCA of batch corrected data all genes.svg')
        else  paste('2D PCA of batch corrected data top 500 genes.svg')
      }
      else if(input$plotType == "3D")
      {
        if(input$top==2)  paste('3D PCA of batch corrected data all genes.pdf')
        else  paste('3D PCA of batch corrected data top 500 genes.pdf')
      }
    },
    
    content = function(file) {
      point_size=1
      req(input$pca_2d_plot_point_size)
      if(input$pca_2d_plot_point_size!="") point_size=as.numeric(input$pca_2d_plot_point_size)
      
      p<-pcaplot(data = data, annotation = input$annotation_pca, top = input$top,
                 plotType = input$plotType, point_size = point_size)[[2]]
      ggsave(file,p, width = 22, height = 15, units = "cm")
    }
  )
  return(reactive({input$top}))
}
