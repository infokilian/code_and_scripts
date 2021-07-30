##When the remove samples button is clicked in the normalized tab
#####Display interactive boxplot
boxplot_output<-function(edata,pData,col)
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  data<-edata
  pheno<-pData


  #Rename sample id with row just for display
  i<-1:nrow(pheno)
  
  #get which column u want to plot by
  var<-as.numeric(col)

  # #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]#1
  colnames(data)<- as.vector(pheno[,var])
  boxplot.matrix(as.matrix(data),outline=FALSE,xlab='Rows',ylab='Value',col=colors[pheno[,var]],boxwex=0.25,names=lapply(i,function(x) paste('row',x[1])))
  legend("topright", legend = unique(pheno[,var]),xpd=TRUE, pch = 16, col = colors[unique(pheno[,var])], cex=0.85,inset=0.0005)

}
   
pca_components<-function(dataset,column_no,annotation,resize_factor=NULL){
  #obtain PC (principle components)
  pca = prcomp(t(dataset),scale=F)
  #signif rounds off the percentage calculated by 3 digits
  pcaVars2=signif(((pca$sdev^2))/(sum((pca$sdev^2))),3)*100
  signed = ifelse(max(pca$x[,2] > 70), 1, -1) # make same sign
  signed1 = ifelse(max(pca$x[,3] > 70), 1, -1) # make same sign
  total_variance = sum(pcaVars2[1:3])
  c<-annotation[,column_no]

  return(list(pca$x[,1],signed*pca$x[,2],signed1*pca$x[,3],pcaVars2[1],pcaVars2[2],pcaVars2[3],c))
  #o/p returned is a list of (PC1,PC2,PC3,Percentage of variance contributed by PC1,
  #                                       Percentage of variance contributed by PC2,
  #                                       Percentage of variance contributed by PC3,
  #                                       annotation)
}

pcaplot<-function(data,annotation,top,plotType,resize_factor=NULL,point_size)
{
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Making plot", value = 0)

  # Number of times we'll go through the loop to update progres bar
  n1 <- 3
  #========user input=======#
  pheno<-data[[2]] #annotation table
  dat<-data[[3]] #rlog transformed data
  var<-as.numeric(annotation)#Get which column u want to plot by
  #========================#

  #compute variance of all genes. this is done by taking the row variance (rowvars function) for each gene.
  rv <- rowVars(dat)
  #pick the top 500 most variable genes. This list will be used if user wants PCA of top 500 genes to be displayed
  genes<-order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]
 
  # Increment the progress bar once input recieved, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 1,"/",n1))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Get the components of PCA
  parameters<-list()
  if(top==2) parameters<-pca_components(dat,var,pheno)# pca components computed for all genes #input
  else parameters<-pca_components(dat[genes,],var,pheno)#pca compnents of top 500 most variable genes 

  # Increment the progress bar after obtaining the components, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 2,"/",n1))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Define color pallete. This is needed to color the sample points on the pca plot.
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  #Obtain the percentage of variance explained by each component
  xlab = paste("PC1:",parameters[[4]],"%")
  ylab = paste("PC2:",parameters[[5]],"%")
  zlab = paste("PC3:",parameters[[6]],"%")
  #condense the principle component vectors(PCA1,PCA2,PCA3) and colors assigned to each sample into a dataframe
  #df is matrix with rows correspond to samples and columns (PCA1,PCA2,PCA3,Colors assigned to each sample)
  df<-data.frame(PCA1=parameters[[1]],PCA2=parameters[[2]],PCA3=parameters[[3]],an=parameters[[7]])
  p<-plotly_empty() %>% layout(autosize = F,width=1000,height=800)
  if(!is.null(resize_factor)) p<-plotly_empty() %>% layout(autosize = F,width=500,height=400)
  #check if the color column of df is not empty
  print("line 105 -functions")
  p0<-NULL
  if(length(colors[unique(df$an)])!=0)
  {

    if (identical(plotType, "2D")) { #variable plotType is input obtained from user: 2d/3d
      print(colors[unique(df$an)])
      p0 <- ggplot(df, aes(x = PCA1, y = PCA2, color=df$an)) + 
        geom_point(size=point_size)+
        labs(x=xlab,y=ylab)+
        scale_color_manual(name=colnames(pheno)[var],values=colors[unique(df$an)])+ theme_bw() #2d PCA plot

      if(!is.null(resize_factor)) p<-ggplotly(p0) %>% layout(autosize = F,dragmode ="select",width=500,height=400)
      else p<-ggplotly(p0) %>% layout(autosize = F,width=1000,height=800,dragmode = "select")
    } else {

      if(!is.null(resize_factor))
      {
        p<- plot_ly(x=df$PCA1,y=df$PCA2,z=df$PCA3) %>%    #3D PCA plot
          add_markers(type = "scatter3d",color=df$an,colors = unique(colors[df$an])) %>%
          layout(autosize = F,width=500,height=400,scene = list(xaxis = list(title = xlab),
                                                                 yaxis = list(title = ylab),
                                                                 zaxis = list(title = zlab)))
      }
      else
      {
        p<- plot_ly(x=df$PCA1,y=df$PCA2,z=df$PCA3) %>%    #3D PCA plot
          add_markers(type = "scatter3d",color=df$an,colors = unique(colors[df$an])) %>%
          layout(autosize = F,width=1000,height=800,scene = list(xaxis = list(title = xlab),
                                                                 yaxis = list(title = ylab),
                                                                 zaxis = list(title = zlab)))
      }
      
    }
  }
  else
    {
      if(!is.null(resize_factor)) p<-plotly_empty() %>% layout(autosize = F,width=500,height=400)
      else p<-plotly_empty() %>% layout(autosize = F,width=1000,height=800)
    }

  # Increment the progress bar when PCA plot has been successfully generated, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 3,"/",n1))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  return(list(p,p0)) #return PCA plot
}

# The following function ensures each treatment/condition has a min of 3 samples 
# else returns treatment/condition groups containing either one or two samples
min_samples_three <- function(dds.norm) {
  print('Inside function min_samples_three')
  dt1<-as.vector(sapply( levels(dds.norm$condition), function(lvl)
  {
    #Obtain number of samples present in each treatmen/conditiont group
    num<-ncol(counts(dds.norm,normalized=TRUE)[,dds.norm$condition == lvl] )
    if(is.null(num)) num<-1
    num
    
  }))
  #dt1 is a dataframe with two columns namely treatment/condition group and number of samples in a treatment group
  #Which treatment/condition groups in dt1 have 2 samples (Obtain index)
  two_samples<-which(dt1 %in% 2)
  #Which treatment/condition groups in dt1 has 1 sample (Obtain index)
  one_sample<-which(dt1 %in% 1)

  #Get a list of all treatment/condition groups
  conditions<-levels(dds.norm$condition)
  
  #return only those treatment/condition groups having only one or two samples
  if(length(one_sample)!=0) return (list(conditions[one_sample],1))#return treatment/condition groups containing only one sample
  else if(length(two_samples)!=0) return (list(conditions[two_samples],2))#return treatment/condition groups containing only two samples
  else if(length(one_sample)!=0 && length(two_samples)!=0) NULL #return null if there are no treatment/condition groups with neither one or two samples
  
  
}

#The following plots (density and ECDF plots) is used as a quality control for normalized data
#It is triggered by clicking the QC button in normalized table tab(represented as a sub tab under normalization tab)
#Plot 1: Density plot function: returns density plot
density_plot<-function(normal,x,y,zoom) #input(normalized table,x and y limits of plots if parameter zoom =True
                                        # Zoom parameter is true when user wants to zoom in,this results in resizing the plot )
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  library(geneplotter)
  #To assess whether the normalization has worked,
  #we plot the densities of counts for the different samples.
  #Since most of the genes are (heavily) affected by the experimental conditions,
  #a succesful normalization will lead to overlapping densities.
  if(zoom==F) multidensity( normal,xlab="mean counts",xlim=c(0,1000),col=colors)
  else multidensity( normal,xlab="mean counts",xlim=x,ylim=y,col=colors)
  
}
#plot 2: ECDF plot function
ECDF_plot<-function(normal)
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  library(geneplotter)
  #In an ECDF plot, the estimated probility is plotted on the y-axis and the count values
  #on the x-axis. I.e. you can read of the median and other quantiles from this plot.
  #As already mentioned, if the normalization has worked, the ECDFs of the different samples
  #should be overlapping.
  multiecdf(normal, normalized = T,xlab="Normalized gene counts", xlim = c(0,1000),col=colors)
  
  }
#batch correction function
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
heatmap_sample<-function(sampleDists,vsd)
{
  library("RColorBrewer")
  library("pheatmap")
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Making plot", value = 0)
  
  # Number of times we'll go through the loop
  n <- 3
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
  colnames(sampleDistMatrix) <- paste(paste(vsd$condition, sep="-"),"(",colnames(vsd),")",sep=" ")#NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
 p<-pheatmap::pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
 # Increment the progress bar, and update the detail text.
 progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
 
 # Pause for 0.1 seconds to simulate a long computation.
 Sys.sleep(0.1)
 dev.off()
 return(p)
}

#ma plot
ma_plot<-function(ma_choice,combination,DE_genes,scale,p_values)
{
  if(!is.null(ma_choice))
  {
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Processing Data", value = 0)
    
    num<- length(combination)
    
    #Get the deseq2 dataset
    res<-DE_genes[[as.numeric(ma_choice)]][5] [[1]]

    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 1,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    #filter most significant genes with fdr cut off 0.05
    resSig <- subset(res, padj < p_values[[as.numeric(ma_choice)]])

    #get down regulated genes
    d<-resSig[resSig$log2FoldChange<1,]
    topGene_down<-rownames(head(d[order(d$log2FoldChange),],2))
    #get up regulated genes
    u<-resSig[resSig$log2FoldChange>1,]
    topGene_up<-rownames(head(u[order(u$log2FoldChange,decreasing = TRUE),],2))
    #top variable genes
    topGene <- rownames(resSig)[which.min(resSig$padj)]# only most significant gene
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 2,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)

    plotMA(res,ylim=c(-1*scale,scale))

      # Increment the progress bar, and update the detail text.
      progress$inc(1/3, detail = paste("Doing part", 3,"/",3))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
  }
}

#volcano plot
volcano_plot<-function(vol_choice,combination,DE_genes,scale_volx,scale_voly,hallmark,GO,organism,dataset,user_gene,user_list,max_num)
{
  # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Processing Data", value = 0)

  #when user specified genes (user_gene) or a user specified gene list (user_list) is selected, annotations will not be filtered 
  #by hallmark pathways or GO terms. Hallmark and Go are always checked if it is available with the choosen organism (is.null(GO), etc.)
   
    if(user_gene != "" & !is.null(user_gene)){ 
      user_gene <- unlist(strsplit(user_gene, ","))
      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]])
      result_up <- result_up[result_up$Row.names %in% user_gene,]
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      result_down <- result_down[result_down$Row.names %in% user_gene,]
    
    } else if(!is.null(user_list)){
      user_list <- as.vector(t(read.delim(file = user_list$datapath, sep = ",", header = F)))
      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]])
      result_up <- result_up[result_up$Row.names %in% user_list,]
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      result_down <- result_down[result_down$Row.names %in% user_list,]
    
    } else if((hallmark != "Default" & GO == "" & !is.null(hallmark)) | (hallmark != "Default" & is.null(GO) & !is.null(hallmark))){

      c1_hallmark <- as.data.frame(msigdbr(species = organism, category = "H"))
      c1_hallmark <- c1_hallmark[c1_hallmark$gs_name==hallmark, "gene_symbol"]
      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]])
      result_up <- result_up[result_up$Row.names %in% c1_hallmark,]
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      result_down <- result_down[result_down$Row.names %in% c1_hallmark,]
    
    } else if((hallmark == "Default" & GO != "" & !is.null(GO)) | (is.null(hallmark) & GO != "" & !is.null(GO))){ 

      mart = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
      mart <- getBM(attributes = c("external_gene_name","go_id"), mart = mart)
      mart <- mart[mart$go_id == GO, "external_gene_name"]
      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]]) 
      result_up <- result_up[result_up$Row.names %in% mart,]
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      result_down <- result_down[result_down$Row.names %in% mart,]
      
    } else if(hallmark != "Default" & GO != "" & !is.null(hallmark) & !is.null(GO)){ 

      mart = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
      mart <- getBM(attributes = c("external_gene_name","go_id"), mart = mart)
      mart <- mart[mart$go_id == GO, "external_gene_name"]
      c1_hallmark <- as.data.frame(msigdbr(species = organism, category = "H"))
      c1_hallmark <- c1_hallmark[c1_hallmark$gs_name==hallmark, "gene_symbol"]
      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]]) 
      result_up <- result_up[result_up$Row.names %in% c(mart,c1_hallmark),]
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      result_down <- result_down[result_down$Row.names %in% c(mart,c1_hallmark),]
      
    } else {

      result_up <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[1]])
      result_down <- as.data.frame(DE_genes[[as.numeric(vol_choice)]][[2]]) 
      
    }
    
  #based on the filtered table of DE genes, labels are chosen based on the specified maximal number of genes (max_num)  
    upDE <-  result_up
    FClabel_up <- upDE[order(abs(upDE$log2FoldChange), decreasing = TRUE),]
    if(nrow(FClabel_up)>max_num){
      FClabel_up <- as.character(FClabel_up[c(1:max_num),"Row.names"])
    } else {
      FClabel_up <- as.character(FClabel_up$Row.names)}
    plabel_up <- upDE[order(upDE$padj, decreasing = FALSE),]
    if(nrow(plabel_up)>max_num){
      plabel_up <- as.character(plabel_up[c(1:max_num),"Row.names"])
    } else {
      plabel_up <- as.character(plabel_up$Row.names)}
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 1,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    downDE <-  result_down
    FClabel_down <- downDE[order(abs(downDE$log2FoldChange), decreasing = TRUE),]
    if(nrow(FClabel_down)>max_num){
      FClabel_down <- as.character(FClabel_down[c(1:max_num),"Row.names"])
    } else {
      FClabel_down <- as.character(FClabel_down$Row.names)}
    plabel_down <- downDE[order(downDE$padj, decreasing = FALSE),]
    if(nrow(plabel_down)>max_num){
      plabel_down <- as.character(plabel_down[c(1:max_num),"Row.names"])
    } else {
      plabel_down <- as.character(plabel_down$Row.names)
    }
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 2,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    #labels are joined in one vector and informations on the labels is added to the data.frame "data"
    label<- unique(c(FClabel_up, plabel_up, FClabel_down, plabel_down))
    
    data <- DE_genes[[as.numeric(vol_choice)]][[3]]
    label_data <- c()
    for (i in 1:nrow(data)){
      if(data$Row.names[i] %in% label){label_data[i] =  as.character(data$Row.names[i])}
      else{ label_data[i] = ""}
    }
    data$label<- label_data

    regulation_data <- c()
    for (i in 1:nrow(data)){
      if(data$log2FoldChange[i]> scale_volx && data$padj[i]< scale_voly){ regulation_data[i] <- "up"}
      else if (data$log2FoldChange[i]< -scale_volx && data$padj[i] < scale_voly){ regulation_data[i] <- "down"}
      else { regulation_data[i] <- "ns"}
    }
    data$regulation = regulation_data
    
    # Volcano Plot - still an issue with the text size of the legend
    g<-ggplot(data=na.omit(data), aes(x=log2FoldChange, y=(-log10(padj)), colour=regulation)) +
      geom_point(alpha=0.4, size=1.75) +
      scale_color_manual(labels=c(paste0(length(regulation_data[regulation_data=="down"])," down-regulated genes"), 
                                paste0(length(regulation_data[regulation_data=="ns"])," non-significant genes"),
                                paste0(length(regulation_data[regulation_data=="up"])," up-regulated genes")), 
                         values=c("cornflowerblue","grey", "firebrick"))+
      scale_x_continuous() +
      scale_y_continuous() +
      xlab("log2(FoldChange)") +
      ylab("-log10(padj)") +
      geom_vline(xintercept = 0, colour="black")+
      geom_vline(xintercept = c(-scale_volx,scale_volx), colour="red")+
      geom_hline(yintercept=-log(scale_voly,10),colour="red")+
      geom_text_repel(data=na.omit(data[!data$label =="",]),aes(label=label), colour = "black", size=4)+
      guides(colour = guide_legend(title = ""))+
      #theme(legend.text = element_text(size = 4))+
      theme_bw()
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 3,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    return(g)
    #g 
}
get_pheno<-function(data,pheno)
{
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
    return(pheno)
  }
}
#Enrichment function(returns object of enrichment function called)
#if Kegg anallysis is performed the enrichKegg function is called and its o/p is returned
enrichment_function<-function(enrichment_type,enrichment_input)
{
  
  obj<-NULL
    if(enrichment_type=="kegg")
    {
      obj<-enrichKEGG(gene         = enrichment_input[[1]],
                      organism     = enrichment_input[[2]],
                      pAdjustMethod = "none",
                      qvalueCutoff = 1,
                      minGSSize = 0)
    }
    else if(enrichment_type=="biological process")
    {
      obj<-enrichGO(gene     = enrichment_input[[1]],
               universe      = enrichment_input[[2]],
               OrgDb         = enrichment_input[[3]],
               ont           = "BP",
               pAdjustMethod = "none",
               qvalueCutoff  = 1, 
               readable      = TRUE)
      
    }
    else if(enrichment_type=="hallmark") 
    {
      
      obj<-enricher(enrichment_input[[1]],
                    TERM2GENE = enrichment_input[[2]],
                    universe  = enrichment_input[[3]],
                    pAdjustMethod = "none",
                    qvalueCutoff = 1,
                    minGSSize = 5)
    }
    
}


enrichment_main<-function(enrichment_type,result,input_organism,dds.fc,num,CoCena,c1_hallmark,orgDB)
{
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Processing data", value = 0)
  
  #get all de genes
  #result
  org<-NULL
  organism<-NULL
  universe<-NULL
  All_genes_symbol<-NULL
  All_genes_ensembltrans<-NULL
  All_genes<-NULL

  #organism information
    org<-orgDB
    
    if(enrichment_type=="kegg") organism<- search_kegg_organism(input_organism, ignore.case = T)[1,1]
    
    else 
    {
     All_genes=AnnotationDbi::select(org,rownames(assay(dds.fc)),"ENTREZID","SYMBOL",multiVals='first')
      if (enrichment_type == "hallmark"){
        universe<- All_genes[!is.na(All_genes$ENTREZID),2]
      }
    }

  n<-num*2
  Enriched_list<-list()
  Enriched_obj<-list()
  Enriched_Kegg_gene<-list()#needed to display user selecetd kegg batch way.(contains a list of genes
                            #pertaining to the selected pathway and the corresponding foldchange values)
  k<-0
  
#if CoCena is performed then we perform enrichment analysis on the 
if(!is.null(CoCena))
    {
      cluster_info<-CoCena
      n<-n+nrow(cluster_info)
    }

for(i in 1:num) # looping through all DE comparisons
  {
    e_list<-list()
    e_obj<-list()
    kegg_genelist<-list() #used only for kegg
    for(j in 1:2) #two loops. First loop for up regulated de genes in a comparison 
    {                         #Second loop for down regulated de genes in a comparison
      k<-k+1
      res<-as.data.frame(result[[i]][j]) #get up/down regulated de genes for a comparison
      if(nrow(res)!=0)
      {
        genes<-res[,1]
        df<-res[,-1]
        rownames(df)<-genes
        entrez_id<-NULL
        gene_symbol<-NULL
        obj<-NULL
        geneList<-NULL #only popuated when kegg is computed

if ((enrichment_type=="kegg") && (!(input_organism %in% c("Homo sapiens", "Mus musculus"))))
          {
          genes_de <- as.character(rownames(df))
          genes_de_uni <- bitr(genes_de, fromType = "SYMBOL", toType = "UNIPROT",org)
          entrez_id <- bitr_kegg(genes_de_uni[,2], fromType = "uniprot", toType = "kegg",organism)[,2]
          entrez_id <- entrez_id[!is.na(entrez_id)]
          reg=AnnotationDbi::select(org,rownames(df),"ENTREZID","SYMBOL",multiVals="first")
          idx <- match(rownames(df), reg$SYMBOL)
          df$entrez<-reg[idx,]$ENTREZID
          gene_symbol<-reg$SYMBOL[idx]
    
          } else {
          #convert de gene names to entrez id(input format required to perform enrichment analysis)
          reg=AnnotationDbi::select(org,rownames(df),"ENTREZID", "SYMBOL",multiVals="first")
          idx <- match(rownames(df), reg$SYMBOL)
          df$entrez<-reg[idx,]$ENTREZID
          gene_symbol<-reg$SYMBOL[idx]
          entrez_id<-df$entrez[!is.na(df$entrez)]

        }
        if(enrichment_type=="kegg")
          {
           input<-list(entrez_id,organism)
           obj<-enrichment_function("kegg",input)
           #compute gene list for kegg
           temp<-data.frame(entrez=c(reg[idx,]$ENTREZID),foldchange=c(df[idx,]$log2FoldChange))

           temp<-temp[complete.cases(temp), ]
           geneList<-temp$foldchange
           names(geneList)<-temp$entrez
          }
        else if(enrichment_type=="hallmark")
          {
            input<-list(entrez_id,c1_hallmark,universe)
            obj<-enrichment_function("hallmark",input)
          }
        else if(enrichment_type=="biological process") 
        {
          input<-list(entrez_id,All_genes$ENTREZID,org)
          obj<-enrichment_function("biological process",input)
        }

        df_obj<-NULL
        if(!is.null(obj))
        {
          if(nrow(as.data.frame(summary(obj)))>0)
          {
            df_obj<-as.data.frame(summary(obj))[1:8]
            if(enrichment_type!="biological process")
            {
              df_obj<-as.data.frame(summary(obj))[1:8]

              for(x in 1:length(df_obj[,8]))
              {

                temp<-strsplit(df_obj[x,8],"/")

                id<-which(entrez_id %in% temp[[1]] )

                df_obj[x,8]<-paste(gene_symbol[id], collapse = '/')
            }
            
            }
            e_list[[length(e_list)+1]]<-df_obj
            e_obj[[length(e_obj)+1]]<-obj
            kegg_genelist[[length(kegg_genelist)+1]]<-geneList
          }
          else
          {
            e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
            e_obj[[length(e_obj)+1]]<-NULL
            kegg_genelist[[length(kegg_genelist)+1]]<-NULL
          }
          
        }
        else
        {
          e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          e_obj[[length(e_obj)+1]]<-NULL
          kegg_genelist[[length(kegg_genelist)+1]]<-NULL
        }
      }
      else
      {
        e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
        e_obj[[length(e_obj)+1]]<-NULL
        kegg_genelist[[length(kegg_genelist)+1]]<-NULL
      }
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", k,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
    e_obj[[length(e_obj)+1]]<-NULL
    kegg_genelist[[length(kegg_genelist)+1]]<-NULL
    Enriched_list[[i]]<-e_list
    Enriched_obj[[i]]<-e_obj
    Enriched_Kegg_gene[[i]]<-kegg_genelist
}
  
if(!is.null(CoCena)) 
    {
      #CoCena
    cluster_info<-CoCena
      for(i in 1:nrow(cluster_info))
      {
        DEG<-unlist(strsplit(cluster_info$gene_n[i], split = ","))
        print(str(DEG))
        e_list<-list()
        e_obj<-list()
        k<-k+1

        if(!identical(character(0),DEG))
        {
          entrez_id<-NULL
          gene_symbol<-NULL
          obj<-NULL
          
          if ((enrichment_type=="kegg") && (!(input_organism %in% c("Homo sapiens", "Mus musculus"))))
          {
            genes_de_uni <- bitr(DEG, fromType = "SYMBOL", toType = "UNIPROT",org)
            entrez_id <- bitr_kegg(genes_de_uni[,2], fromType = "uniprot", toType = "kegg",organism)[,2]
            entrez_id <- entrez_id[!is.na(entrez_id)]
            reg=AnnotationDbi::select(org,DEG,"ENTREZID","SYMBOL",multiVals="first")
            idx <- match(DEG, reg$SYMBOL)
            entrez<-reg[idx,]$ENTREZID
            gene_symbol<-reg$SYMBOL[idx]
            
          } else {
            #convert de gene names to entrez id(input format required to perform enrichment analysis)
            reg=AnnotationDbi::select(org,DEG,"ENTREZID", "SYMBOL",multiVals="first")
            idx <- match(DEG, reg$SYMBOL)
            entrez<-reg[idx,]$ENTREZID
            gene_symbol<-reg$SYMBOL[idx]
            entrez_id<-entrez[!is.na(entrez)]
            
          }
          if(enrichment_type=="kegg")
          {
            input<-list(entrez_id,organism)
            obj<-enrichment_function("kegg",input)

          }
          else if(enrichment_type=="hallmark")
          {
            input<-list(entrez_id,c1_hallmark,universe)
            obj<-enrichment_function("hallmark",input)
          }
          else if(enrichment_type=="biological process") 
          {
            input<-list(entrez_id,All_genes$ENTREZID,org)
            obj<-enrichment_function("biological process",input)
          
          }
          df_obj<-NULL
          if(!is.null(obj))
          {
            if(nrow(as.data.frame(summary(obj)))>0)
            {
              df_obj<-as.data.frame(summary(obj))[1:8]
              if(enrichment_type!="biological process")
              {
                for(x in 1:length(df_obj[,8]))
                {
                  temp<-strsplit(df_obj[x,8],"/")
                  id<-which(entrez_id %in% temp[[1]] ) #get entrez id of genes in a pathway
                  df_obj[x,8]<-paste(gene_symbol[id], collapse = '/')#convert entrez id to genes

                }
              }
              
              d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
              Enriched_list[[length(Enriched_list)+1]]<-list(d,d,df_obj)
              Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,obj)
              Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
            }
            else
            {
              d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
              Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
              Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
              Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
            }
            
          }
          else
          {
            d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
            Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
            Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
            Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
          }
        }
        else
        {
          d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
          Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
          Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
        }
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", k,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    }
  
  return(list(Enriched_list,Enriched_obj,Enriched_Kegg_gene))
}

#Display barplot of top 10 kegg pathway identified for selected comparison

enrichment_plot<-function(enrichment_type,result,res,row,col,category,category_go){
  
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
  
  #up regulated column
  up<-unique(nrow(res[[row]][[1]]))

  #down regulated column
  down<-unique(nrow(res[[row]][[2]]))
  
  #which column is null
  col_idx<-NULL
  if(up==0) {col_idx<- 1}
  else if (down==0) {col_idx<-2}

  
  if(enrichment_type!="biological process"){

      kk<-result[[row]][[col]]
      p<-barplot(kk, showCategory=as.numeric(category))
    
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      return(p)

  }
  
  else{
    if(nrow(res[[row]][[col]])!=0) ego<-result[[row]][[col]]
    else ego<-NULL 
    
    x<-NULL
    if(category_go!=""){
      x<-gofilter(ego,level=as.numeric(category_go))#dropGO(ego, level = as.numeric(input$category_go), term = NULL)
    }
    else x<-ego
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    if(!is.null(x)){
        p<-barplot(x, showCategory=as.numeric(category))
    }
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    return(p)
  }
}
####Compute transcription factors from a list
Enriched_transcription_factors<-function(TF_list,result,num,input_organism,CoCena,
                                         anova_table,conchoice,dds.fc)
{
  if (!is.null(TF_list)){  
    #get all de genes
    
    DE_TF<-list()
    for(i in 1:num)
    {
      TF<-list()
      for(j in 1:3)
      {
        res<-as.data.frame(result[[i]][j])
        genes<-res[,1]
        df<-res[,-1]
        rownames(df)<-genes

          TF[[length(TF)+1]]<-df[which(rownames(df) %in% TF_list),]

      }
      DE_TF[[i]]<-TF
    }
    if(!is.null(CoCena))
    {
      cluster_info<-CoCena
      #########preparing the anova table in the order as output#########
      a_tab<-anova_table[,-c(2,3)]
      cond<-unique(colData(dds.fc)[,as.numeric(conchoice)])
      c<-colnames(a_tab)
      temp<-as.vector(c[4:(3+length(cond))])
      temp2<-as.vector(c[(length(cond)+4):length(c)])
      all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
      ##################################################################
      for(i in 1:nrow(cluster_info))
      {
        
        DEG<-unlist(strsplit(cluster_info$gene_n[i], split=","))
          genelist<-DEG[which(DEG %in% TF_list)]
          anova_genes<-rownames(all_genes)
          df_final<-all_genes[which(anova_genes %in% genelist),]
          DE_TF[[num+i]]<-list(NULL,NULL,df_final)
        
      }
      
    }
    return(DE_TF)
  }
}
#Heatmap function for either DE genes/DE transcription factors and ANOVA
heatmap_genes<-function(heatmap_call,dds,dds.fc,rld,heat_choice,
                        DE_genes,cluster_info,num,heatmap_name,
                        Distance,Linkage)
{
  
  topVarGenes=NULL
  print("function line 1186")

  if(heatmap_call=="ANOVA")
  {
    resSig<-results(dds)
    
    topVar<-head(rownames(resSig[order(resSig$padj),]),1000)
    topVarGenes<-which(rownames(rld) %in% topVar)
    
  }
  else if(heatmap_call=="DE")
  {
    if(!is.null(heat_choice))
    {
      result<-DE_genes
      DEG<-NULL
      
     if(as.numeric(heat_choice)>num)
      {
        DEG<-unlist(strsplit(cluster_info$gene_n[as.numeric(heat_choice)-num], split=","))
      }

        DEG<-as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]
        
      topVarGenes<-which(rownames(rld) %in% DEG)

    }
    
  }
  else if(heatmap_call=="TF")
  {
    if(!is.null(heat_choice))
    {
      result<-DE_genes

      DEG<-NULL
      
      if(as.numeric(heat_choice)>num)
      {
        
        DEG<-result[[as.numeric(heat_choice)]][[3]]
      }
      
      

        DEG<-rownames(result[[as.numeric(heat_choice)]][[3]])
      

      topVarGenes<-which(rownames(rld) %in% DEG)

    }
    
  }else if(heatmap_call=="TF_custom")
  {
    if(!is.null(heat_choice))
    {
      result<-DE_genes
      DEG<-NULL
      
      if(as.numeric(heat_choice)>num)
      {
        
        DEG<-result[[as.numeric(heat_choice)]][[3]]

      }
      else
      {

        DEG<-rownames(result[[as.numeric(heat_choice)]][[3]])
      }
      topVarGenes<-which(rownames(rld) %in% DEG)

    }
    
  }
  shiny::validate(need(!is.null(topVarGenes),"The number of differentially expressed genes for this condition is Zero"))
  library(RColorBrewer)
  my_palette <- colorRampPalette(c("blue","white", "red"))(n = 15)

  #Color key
  library(RColorBrewer)
  colors <-brewer.pal(8, "Dark2")
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]

  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  mat<-rld[topVarGenes,]
  p<-NULL
  if(nrow(mat)>1 && !is.null(nrow(mat)))
  {

  mat <- t(scale(t(mat))) # scale and center rows
  
  dist_method<-c("euclidean", "manhattan")

  distance = dist(mat,method = dist_method[as.numeric(Distance)])
  link_method<-c("average", "complete","Ward.D2","Ward.D","single")
  cluster = hclust(distance, method = link_method[as.numeric(Linkage)])
  

  names<-FALSE
  if(nrow(mat)<20) names<-TRUE

  annotation <- data.frame(condition = colData(dds.fc)[,"condition"])
  rownames(annotation) <-  colnames(mat)# check out the row names of annotation
  
  Var1        <- unique(colors[colData(dds.fc)[,"condition"]])
  names(Var1) <- levels(colData(dds.fc)[,"condition"])
  anno_colors <- list(condition = Var1)

  p<-pheatmap::pheatmap(mat, annotation = annotation,color = my_palette,
              main = heatmap_name,
              clustering_distance_rows = dist_method[as.numeric(Distance)],
              clustering_method = link_method[as.numeric(Linkage)],
              annotation_colors = anno_colors,
              show_rownames=names)
}
  return(p)
  }

  
#function to check for empty lists
has_empty_list <- function(x) {
  if(is.list(x)) {
    if (length(x)==0) {
      return(TRUE)
    } else {
      return(any(vapply(x, has_empty_list, logical(1))))
    }
  } else {
    return(FALSE)
  }
}

#source of variation output
source_of_variation_op<-function(edata,pData,ids)
{
  sov<-as.data.frame(matrix(NA,nrow=length(ids),ncol=2))
  colnames(sov)<-c("Variation","Percentage")
  for(i in 1:length(ids))
  {
    sov[i,1]<-colnames(pData)[ids[i]]
    sov[i,2]<-source_of_variation_ip(edata,pData,ids[i])
  }

  return(sov)
} 

source_of_variation_ip<-function(edata,pData,k)
{
  y_bar<-sum(rowSums(edata))/length(edata)

  #SSB(Sum of square between groups) and
  #SSW(Sum of squares within each group) calculation
  SSB=0
  SSW=0
  temp<-c()
  print(levels(pData[,k]))
  #k is the annotation variable
  #Loop through each group in k
  for(i in levels(pData[,k]))
  {
    #get which samples belong to a particular group
    id<-which(pData[,k]==i)
    #group mean
    y_i<-mean((rowSums(edata[,id])))
    group<-0
    
    #loop through each sample in a group
    for(j in id)
    {
      p<-sapply(edata[,j],function(idx) (y_i-idx)^2)
      group<-group+sum(p)

    }
    SSW<-SSW+group
    temp<-c(temp,(y_i-y_bar)^2)
  }
  SSB<-sum(temp)*length(edata)
  TSS<-SSB+SSW
  R<-SSB/TSS
  R_Square=R^2

  return(R_Square)
}

plotCounts<-function(batch_corrected=NULL,dds, gene, intgroup = "condition", normalized = TRUE, 
                     transform = FALSE, main, xlab = "group", returnData = FALSE, 
                     replaced = FALSE)
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
                                                             "factor"))))
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene,]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  data <- data.frame(count = cnts + 0.5, group = as.integer(group))
  
  if(!is.null(batch_corrected))
  {
    cnts <- batch_corrected
    data <- data.frame(count = cnts, group = as.integer(group))
  }
  if (transform) {
    data$count <- log2(data$count)
    ylab <- expression(log[2] ~ count)
    logxy <- ""
  }
  else {
    ylab <- ifelse(normalized, "normalized count", "count")
    logxy <- "y"
  }
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  if (returnData)
    return(data.frame(count = data$count, colData(dds)[intgroup]))
  plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count,
       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
       xlab = xlab, ylab = ylab, main = main)
  axis(1, at = seq_along(levels(group)), levels(group))
}

gene_counts<-function(dataset,dds.fc,gene_name,log_scale)
{
  geneCounts <- plotCounts(dataset,dds.fc, gene=gene_name, intgroup="condition", returnData=TRUE)
  comparisons <- combn(levels(geneCounts$condition), 2, simplify = F)
  

  p<-ggplot(geneCounts, aes(x=condition, y=count,fill=condition)) +
    geom_boxplot()+
    geom_point(position=position_jitter(width=.1,height=0), size=3)+
    stat_compare_means(comparisons = comparisons,label = "p.signif",
                       method = "t.test")+
    ggtitle(paste0("Gene name: ",gene_name)) +
    labs(x="Condition",y="Batch corrected counts") +
    theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=30, hjust=0.5)) +
    theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22))+
    theme_bw()
  if(log_scale=="Yes") p+scale_y_continuous(trans='log2')
  else p
}

summary_analysis<-function(result_de,result_tf,result_kegg,result_bp,result_hall,combination)
{

    res_de<-data.frame(matrix(NA, nrow = length(combination), ncol = 3))
    rownames(res_de)<-lapply(1:length(combination), function(i) {
      combination[[i]]
    })
    res_tf<-data.frame(matrix(NA, nrow = length(combination), ncol = 3))
    rownames(res_tf)<-lapply(1:length(combination), function(i) {
      combination[[i]]
    })
    res_bp<-data.frame(matrix(NA, nrow = length(combination), ncol = 2))
    rownames(res_bp)<-lapply(1:length(combination), function(i) {
      combination[[i]]
      
    })
    res_kegg<-data.frame(matrix(NA, nrow = length(combination), ncol = 2))
    rownames(res_kegg)<-lapply(1:length(combination), function(i) {
      combination[[i]]
      
    })
    res_hall<-data.frame(matrix(NA, nrow = length(combination), ncol = 2))
    rownames(res_hall)<-lapply(1:length(combination), function(i) {
      combination[[i]]
      
    })
    colnames(res_de)<-c('up-regulated_genes','down-regulated_genes','both')
    colnames(res_tf)<-c('up-regulated_TF','down-regulated_TF','All DE TF')
    colnames(res_bp)<-c('Enriched BP for Up-reg genes','Enriched BP for Down-reg genes')
    colnames(res_kegg)<-c('Enriched Kegg pathways for Up-reg genes','Enriched Kegg pathways for Down-reg genes')
    colnames(res_hall)<-c('Enriched Hallmark gene sets for Up-reg genes','Enriched Hallmark gene sets for Down-reg genes')

    for(i in 1:length(combination))
    {
      #de
      res_de[i,1]<-nrow(as.data.frame(result_de[[i]][1]))
      res_de[i,2]<-nrow(as.data.frame(result_de[[i]][2]))
      res_de[i,3]<-nrow(as.data.frame(result_de[[i]][3]))
      #tf
      res_tf[i,1]<-nrow(as.data.frame(result_tf[[i]][[1]]))
      res_tf[i,2]<-nrow(as.data.frame(result_tf[[i]][[2]]))
      res_tf[i,3]<-nrow(as.data.frame(result_tf[[i]][[3]]))

      res_bp[i,1]<-nrow(as.data.frame(result_bp[[i]][[1]]))
      res_bp[i,2]<-nrow(as.data.frame(result_bp[[i]][[2]]))

      res_kegg[i,1]<-nrow(as.data.frame(result_kegg[[i]][[1]]))
      res_kegg[i,2]<-nrow(as.data.frame(result_kegg[[i]][[2]]))

      res_hall[i,1]<-nrow(as.data.frame(result_hall[[i]][[1]]))
      res_hall[i,2]<-nrow(as.data.frame(result_hall[[i]][[2]]))
    }
    return(list(res_de,res_tf,res_kegg,res_bp,res_hall))
  }

#all p value plot
p_value_all<-function(result,p_choice,combination)
{
  
  
  res<-as.data.frame(result[[p_choice]][4])

  p<-hist(res$pvalue[res$baseMean > 1],main = paste("Corrected p-value",combination[[p_choice]]),  col = "lavender", xlab = "p-values")
  return(p)
} 


#TF prediction using ChEA3
TF_prediction_ChEA3<-function(DE_genes,dds.fc,anova_table,
                                combination,conchoice,CoCena,
                                organism,dataset, topTF){
# Create a Progress object
progress <- shiny::Progress$new()
# Make sure it closes when we exit this reactive, even if there's an error
on.exit(progress$close())

#get all de genes
result<-DE_genes
combo<-combination

progress$set(message = "Processing Data", value = 0)

Enriched_dt<-list()
All_pred_tf<-list()
All_idx<-list()

#preparing the anova table in the order as output
a_tab<-anova_table[,-c(2,3)]
cond<-unique(colData(dds.fc[[1]])[,as.numeric(conchoice)])

c<-colnames(a_tab)

temp<-as.vector(c[4:(3+length(cond))])
temp2<-as.vector(c[(length(cond)+4):length(c)])

all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]


for(i in 1:length(combo))
{
  dt<-list()
  All_pred<-list()
  idx_pred<-list()

  for(j in 1:3)
  {
    print(j)
    res<-as.data.frame(result[[i]][[j]])
    
    TFs<-NULL
    TF_table<-NULL
    if(nrow(res)!=0)
    {
      genes<-res[,1]
      df<-res[,-1]
      rownames(df)<-genes  
      
      #ChEA3 can automatically use human and mouse genes, other organisms have to be transformed
        #mouse genes have to be changed back to mouse symbols, thus, requiring an exception in the code
      if(organism == "Homo sapiens"){
        url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
        encode = "json"
        payload = list(query_name = "myQuery", gene_set = genes)
        
        #POST to ChEA3 server
        response = POST(url = url, body = payload, encode = encode)
        json = httr::content(response, as = "text")
        
        #results as list of R dataframes
        TF_table = fromJSON(json)
        TF_table <- TF_table$`Integrated--meanRank`
        
        # extract those from meanRank since meanRank scored as best method:
        TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
        TF_table$Score <- as.numeric(TF_table$Score)
        if (topTF < nrow(TF_table)) TF_table = TF_table[1:topTF,]
        TFs <- as.character(TF_table$TF)
        
        
      } else if(organism == "Mus musculus"){
        url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
        encode = "json"
        payload = list(query_name = "myQuery", gene_set = genes)
        
        #POST to ChEA3 server
        response = POST(url = url, body = payload, encode = encode)
        json = httr::content(response, as = "text")
        
        #results as list of R dataframes
        TF_table = fromJSON(json)
        TF_table <- TF_table$`Integrated--meanRank`
        
        # extract those from meanRank since meanRank scored as best method:
        TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
        TF_table$Score <- as.numeric(TF_table$Score)

        if (topTF< nrow(TF_table)) TF_table = TF_table[1:topTF,]
        
        #change back to mouse symobls
        human = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
        mouse = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
        
        mouse_gene = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = TF_table$TF , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
        TFs <- mouse_gene[, 2]
        
        
        
      } else {
        createAlert(session, "Anchor1", "Warning", title = "Warning!",
                    content = "For the chosen organism human transcription factors are returned")
        #change organism specific gene symbols to human gene symbols
        human = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
        chosen_organism = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
        
        human_gene <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = chosen_organism, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        genes <- human_gene[, 2]
        
        url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
        encode = "json"
        payload = list(query_name = "myQuery", gene_set = genes)
        
        #POST to ChEA3 server
        response = POST(url = url, body = payload, encode = encode)
        json = httr::content(response, as = "text")
        
        #results as list of R dataframes
        TF_table = fromJSON(json)
        TF_table <- TF_table$`Integrated--meanRank`
        
        # extract those from meanRank since meanRank scored as best method:
        TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
        TF_table$Score <- as.numeric(TF_table$Score)
        
        if (topTF < nrow(TF_table)) TF_table = TF_table[1:topTF,]
        TFs <- TF_table$TF
        
        #change human gene symbols back to organism specific gene symbols
        chosen_organism_gene = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = TFs , mart = human, attributesL = c("external_gene_name"), martL = chosen_organism, uniqueRows=T)
        TFs <- chosen_organism_gene[, 2]
        
      }
      if(!is.null(TFs))
      {

        dt[[length(dt)+1]]<-all_genes[as.character(rownames(all_genes)) %in% TFs,]
        All_pred[[length(All_pred)+1]]<-TF_table
        de_genes_pred_idx<- rownames(all_genes)[as.character(rownames(all_genes)) %in% rownames(res)]
        idx_pred[[length(idx_pred)+1]]<-de_genes_pred_idx
        
        
      }
      else
      {

        dt[[length(dt)+1]]<-NULL
        All_pred[[length(All_pred)+1]]<-NULL
        idx_pred[[length(idx_pred)+1]]<-NULL
      }
    }
    else
    {

      dt[[length(dt)+1]]<-NULL 
      All_pred[[length(All_pred)+1]]<-NULL
      idx_pred[[length(idx_pred)+1]]<-NULL
    }
  }
  

  Enriched_dt[[i]]<-dt
  All_pred_tf[[i]]<-All_pred
  All_idx[[i]]<-idx_pred

}
if(!is.null(CoCena))
{
  cluster_info<-CoCena
    for(i in 1:nrow(cluster_info))
    {
      DEG<-unlist(strsplit(cluster_info$gene_n[i], split = ","))
      if(!identical(character(0),DEG))
      {
        if(organism == "Homo sapiens"){
          url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
          encode = "json"
          payload = list(query_name = "myQuery", gene_set = DEG)
          
          #POST to ChEA3 server
          response = POST(url = url, body = payload, encode = encode)
          json = httr::content(response, as = "text")
          
          #results as list of R dataframes
          TF_table = fromJSON(json)
          TF_table <- TF_table$`Integrated--meanRank`
          
          # extract those from meanRank since meanRank scored as best method:
          TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
          TF_table$Score <- as.numeric(TF_table$Score)
          if (topTF < nrow(TF_table)) TF_table = TF_table[1:topTF,]
          TFs <- as.character(TF_table$TF)
          
          
        } else if(organism == "Mus musculus"){
          url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
          encode = "json"
          payload = list(query_name = "myQuery", gene_set = DEG)
          
          #POST to ChEA3 server
          response = POST(url = url, body = payload, encode = encode)
          json = httr::content(response, as = "text")
          
          #results as list of R dataframes
          TF_table = fromJSON(json)
          TF_table <- TF_table$`Integrated--meanRank`
          
          # extract those from meanRank since meanRank scored as best method:
          TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
          TF_table$Score <- as.numeric(TF_table$Score)
          
          if (topTF< nrow(TF_table)) TF_table = TF_table[1:topTF,]
          
          #change back to mouse symobls
          human = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
          mouse = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
          
          mouse_gene = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = TF_table$TF , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
          TFs <- mouse_gene[, 2]
          
          
          
        } else {
          human = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
          chosen_organism = useMart(host = "https://nov2020.archive.ensembl.org","ENSEMBL_MART_ENSEMBL", dataset = dataset)
          
          human_gene <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = DEG , mart = chosen_organism, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
          DEG <- human_gene[, 2]
          
          url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
          encode = "json"
          payload = list(query_name = "myQuery", gene_set = DEG)
          
          response = POST(url = url, body = payload, encode = encode)
          json = httr::content(response, as = "text")
          
          #results as list of R dataframes
          TF_table = fromJSON(json)
          TF_table <- TF_table$`Integrated--meanRank`
          
          # extract those from meanRank since meanRank scored as best method:
          TF_table <- TF_table[,c("TF", "Score", "Overlapping_Genes")]
          TF_table$Score <- as.numeric(TF_table$Score)
          
          if (topTF< nrow(TF_table)) TF_table = TF_table[1:topTF,]
          
          chosen_organism_gene = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = TFs , mart = human, attributesL = c("external_gene_name"), martL = chosen_organism, uniqueRows=T)
          TFs <- chosen_organism_gene[, 2]
        } 
        
        if(!is.null(TFs))
        {
          
          dt[[length(dt)+1]]<-all_genes[as.character(rownames(all_genes)) %in% TFs,]
          All_pred[[length(All_pred)+1]]<-TF_table
          de_genes_pred_idx<- rownames(all_genes)[as.character(rownames(all_genes)) %in% rownames(res)]
          idx_pred[[length(idx_pred)+1]]<-de_genes_pred_idx
          
        }
        else
        {
          dt[[length(dt)+1]]<-NULL
          All_pred[[length(All_pred)+1]]<-NULL
          idx_pred[[length(idx_pred)+1]]<-NULL
        }
        
        
      }
      else
      {
        dt[[length(dt)+1]]<-NULL
        All_pred[[length(All_pred)+1]]<-NULL
        idx_pred[[length(idx_pred)+1]]<-NULL
      }
      
      Enriched_dt[[length(Enriched_dt)+1]]<-dt
      All_pred_tf[[length(All_pred_tf)+1]]<-All_pred
      All_idx[[length(All_idx)+1]]<-idx_pred  
    }
}

return(list(Enriched_dt,All_pred_tf,All_idx))
}  

#########CoCena specific functions#######
correlation_actions <- function(dd2, range_corr, min_corr){ 
  
  output <- list()
  #matrix operation valid for small dataset
  #the result is a list of 3 matrices, (r) coefficient, (n) number of samples and (P) values
  #it is not a data frame but a list of 3 df/matrices
  
  correlation_matrix=Hmisc::rcorr(as.matrix(dd2))
  
  #calculate number of combinations possible in pairwise gene vs gene comparisons (n*n-1/2)
  
  #do the bonferroni correction for the matrix p vals (adj pval = pval *number of comparisons)
  #correlation_matrix$P= correlation_matrix$P*combinations
  #create a 4 column data frame node1(gene1),node2(gene2), rval, pval(adj)
  
  correlation_df = as.data.frame(t(combinat::combn(row.names(correlation_matrix$r), m=2)))
  correlation_df$rval = correlation_matrix$r[lower.tri(correlation_matrix$r)]
  correlation_df$pval = correlation_matrix$P[lower.tri(correlation_matrix$P)]
  
  #retain rows which have pval (adj) < 0.05, and correlations above 0
  correlation_df_filt = correlation_df[correlation_df$pval<0.05 & correlation_df$rval>0,]
  
  #range of cutoff min to max (correlation)
  range_cutoff<-seq(from = min_corr , to = max(correlation_df$rval) , 
                    length.out = range_corr)
  range_cutoff<-round(range_cutoff, 3)
  if(length(range_cutoff)> as.numeric(range_corr)) {
    range_cutoff=range_cutoff[1:range_corr]
  } else {
    range_cutoff=range_cutoff
  }
  
  output[["correlation_df"]] <- correlation_df
  output[["correlation_df_filt"]] <- correlation_df_filt
  output[["range_cutoff"]] <- range_cutoff
  
  return(output)
}  
  
#prepare cutoff calculation
cutoff_prep=function(cutoff, corrdf_r, min_nodes_network){

  ##define function for calculating graph statistics
  ##using only slightly modified version of the original algo

  rsquaredfun = function(graph_df, cutoff){

    #igraph object
    filt_igraph = igraph::graph_from_data_frame(graph_df,directed = F, vertices = NULL)
    #number of nodes
    num_nodes=igraph::vcount(filt_igraph)
    #number of edges
    num_edges=igraph::gsize(filt_igraph)
    #number of components
    graph_components <- igraph::components(filt_igraph)
    num_networks = graph_components$csize[graph_components$csize >= min_nodes_network] %>% 
      length()
    
    #calculate stats
    d = igraph::degree(filt_igraph, mode = "all")
    dd = degree.distribution(filt_igraph, mode = "all", cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    
    
    
    if(length(probability)==0){
      R.square<-0
    }else{
      reg = lm(log(probability) ~ log(degree))
      cozf = coef(reg)
      power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
      alpha = -cozf[[2]]
      R.square = summary(reg)$r.squared
    }

    output =  data.frame(R.squared=R.square,
                         degree=degree,
                         Probs=probability,
                         cutoff=cutoff,
                         no_edges=num_edges,
                         no_nodes=num_nodes,
                         no_of_networks=num_networks)
    
    return(output)
    
  }

  
  ###filter correlations above the cutoff
  filteredmatrix = corrdf_r[corrdf_r$rval>cutoff,]

  ##create expected df with initial values
  output=data.frame(R.squared=0,
                    degree=0,
                    Probs=0,
                    cutoff=cutoff,
                    no_edges=0,
                    no_nodes=0,
                    no_of_networks=0)
  
  #create a switch, if number of rows =0 then just output the initial output else do calculations
  #@shobhit check if the length of the nrow zero df and gtzero df are required to be of the same length
  rownums = ifelse(nrow(filteredmatrix)==0, "zero_rows", "gtzero")
  
  switch(rownums, zero_rows ={output},
         gtzero={rsquaredfun(graph_df=filteredmatrix, cutoff = cutoff)})
  
}

#calculate optimal minimal correlation cutoff
optcut_fun <- function(cutoff_stats){
  
  output <- list()
  
  cutoff_stats_concise = cutoff_stats %>% 
    dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>% 
    dplyr::distinct()
  
  cutoff_stats_concise = cutoff_stats_concise %>% filter(no_of_networks!=0)
  rownames(cutoff_stats_concise) = cutoff_stats_concise$cutoff
  cutoff_stats_concise = cutoff_stats_concise[,-2]
  
  
  
  crit_minmax = c("max","max", "max", "min" )
  names(crit_minmax) = colnames(cutoff_stats_concise)

  normalizationTypes <- rep("percentageOfMax", ncol(cutoff_stats_concise))
  names(normalizationTypes) = colnames(cutoff_stats_concise)
  if(nrow(cutoff_stats_concise)==0){
    return(data.frame())
  }
  nPT = normalizePerformanceTable(cutoff_stats_concise[,c("R.squared", "no_edges", "no_nodes", "no_of_networks")], normalizationTypes)
  w = c(0.5,0.1,0.5, -1)
  names(w) <- colnames(nPT)
  ws<-weightedSum(nPT,w)
  ranked_ws <- rank(-ws) %>% sort()

  calculated_optimal_cutoff <- as.numeric(names(ranked_ws[1]))

  stats_calculated_optimal_cutoff <- cutoff_stats[cutoff_stats$cutoff == calculated_optimal_cutoff, c("degree", "Probs")]
 
  stats <- cutoff_stats_concise
  stats$cut_off <- rownames(stats)
  stats <- stats[stats$cut_off == calculated_optimal_cutoff,]
   
  dd_plot_calculated_optimal = ggplot(stats_calculated_optimal_cutoff,aes(x=log(degree), y= log(Probs))) +
    geom_point() +
    geom_smooth(method="lm") +
    theme_bw() + 
    ggtitle(paste0("Calculated optimal correlation cut-off :",calculated_optimal_cutoff,"; R.squared: ", round(stats[1],3), "; no. edges: ",
                   stats[2], "; no. nodes: ", stats[3], "; no. networks: ", stats[4]))
  
  
  output[["cutoff_stats_concise"]] <- cutoff_stats_concise
  output[["dd_plot_calculated_optimal"]] <- dd_plot_calculated_optimal
  output[["optimal_cutoff"]] <- calculated_optimal_cutoff
  return(output)
  
}

#calculate plot of correlation statistics
plot_cutoffs_internal <- function(cutoff_stats, 
                                  hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL))
  {
  cutoff_stats$corr <- rownames(cutoff_stats) %>% as.numeric()
  
  
  p1 <- plot_ly(cutoff_stats, x = ~corr, y = ~R.squared, type = 'scatter', 
                mode = 'lines+markers', name = "R.squared", line = list(color = "lightblue"), marker = list(color = "lightblue")) 
  p2 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_edges, type = 'scatter', 
                mode = 'lines+markers', name = "no. edges", line = list(color = "orange"), marker = list(color = "orange"))
  p3 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_nodes, type = 'scatter', 
                mode = 'lines+markers', name = "no. nodes", line = list(color = "lightgreen"), marker = list(color = "lightgreen"))
  p4 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_of_networks, type = 'scatter', 
                mode = "markers", name = "no. networks", marker = list(color = "yellow"))
  p <- plotly::subplot(p1, p2, p3, p4, nrows = 4, shareX = T)
  
  steps <- list()
  for(i in 1:length(cutoff_stats$corr)){
    
    step <- list(args = list("marker.color",list(rep("lightblue", length(cutoff_stats$corr)),
                                                 rep("orange", length(cutoff_stats$corr)),
                                                 rep("lightgreen", length(cutoff_stats$corr)),
                                                 rep("yellow", length(cutoff_stats$corr)))), 
                 label = paste0(as.character(cutoff_stats$corr[i]), ", R.squared: ", 
                                round(cutoff_stats$R.squared[i], 3),
                                "; no. edges: ", cutoff_stats$no_edges[i], "; no. nodes: ", 
                                cutoff_stats$no_nodes[i], "; no. networks: ", cutoff_stats$no_of_networks[i]), 
                 method = "restyle"
                 
    )
    for(j in 1:4){
      step[["args"]][[2]][[j]][i] <- "red"
    }
    
    steps[[i]] <- step
  }
  
  
  p <- p %>% layout(hovermode = "x unified")%>%
    layout(title = paste0("Cut-off selection guide"),
           sliders = list(
             list(pad = list(t=60),
                  active = 2, 
                  currentvalue = list(prefix = "Cut-off: ", font = list(color = "black", size = 14)), 
                  steps = steps,
                  font = list(color = "white", size = 0))),
           shapes = list(line))
  return(p)
}

###calculate group fold changes
gfc_calc<-function(grp, trans_norm, group_means, range_GFC){
  df1= trans_norm[,grp]
  df2=gtools::foldchange(df1, group_means) %>% ifelse(.>range_GFC, range_GFC,.) %>% 
    ifelse(.< (-range_GFC), -range_GFC,.) %>%
    as.data.frame()
  colnames(df2) = paste0("",grp)
  return(df2)
}

###compute the different modules
cluster_calculations =function(graph_obj, algo, case, iter) {

  cfg= get(algo)(graph_obj)
  
  mod_score =modularity(graph_obj, cfg$membership)
  
  mod_df= data.frame(modularity_score=mod_score, cluster_algorithm=algo, stringsAsFactors = F)
  
  ##making switch so that in the end when only the best algorithm is to be used then the same function can be used
  output = switch(case, best= cfg$membership, test= mod_df, final=cfg)

  return(output)
}

###calculate mean group fold change for each cluster
gfc_mean_clustergene=function(rownum, cluster_df, gfc_dat){
  d1 = cluster_df[rownum,]
  gene_names= d1["gene_n"] %>%
    stri_split_regex(pattern = ",") %>%
    unlist()
  gfc_means = gfc_dat[gfc_dat$Gene%in%gene_names,] %>%
    dplyr::select(-Gene) %>%
    colMeans()
  d1$conditions = paste0(names(gfc_means), collapse = "#")
  d1$grp_means = paste0(round(gfc_means,3) , collapse = ",")
  return(d1)
}
