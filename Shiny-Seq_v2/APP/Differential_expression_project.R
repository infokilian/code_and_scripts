source("gene_count_module.R") 
Module_Differential_Expression_UI<-function(id)
{
  ns<-NS(id)
  tagList(br(),
          fluidRow(column(4,
                            uiOutput(ns("group1")),
                            uiOutput(ns("group2")),
                            actionButton(ns("combo"), label = "Add combination"))),
          fluidRow(column(4,uiOutput(ns("group3")))),
          fluidRow(column(2, uiOutput(ns("comb"))),
                   column(2,uiOutput(ns("sliders"))),
                   column(2, uiOutput(ns("sliders_fc"))),
                   column(2,uiOutput(ns("ind"))),
                   column(4,uiOutput(ns("choices")))),
          fluidRow(column(2, selectInput(ns("shrinkornot")  ,label = h5("Fold Change Shrinkage?"),
                                         choices = list("yes" = 1, "no" = 2),
                                         selected = 1))),
          br(),
          br(),
          hr(),
    actionButton(ns("ok3"), label = "Start", icon = icon("paper-plane")),
    bscols(widths = 10,div(style="height: 15px;",width = '200px')),
    conditionalPanel(sprintf("input['%s']>0",ns("ok3")),
                     actionButton(ns("p-value"),"Display p-value plot")),
    radioButtons(ns("lfc_de"), "Show log2 foldchange?:", choices = c("Yes", "No"),selected = "No"),
    DT::dataTableOutput(ns("de_genes")),
    br(),
    br(),
    hr(),
    fluidRow(
      column(1,
             selectInput(ns("datachoice4")  ,label = h5("Select Data Type"),
                         choices = list("Excel" = 1, "CSV" = 2),
                         selected = 1)),
      column(1,
             br(),
             br(),
             downloadButton(ns('download_DE_genes_Table'), 'Download DE Genes Table'))),
    DT::dataTableOutput(ns("filtered_data")),
    gene_count_module_UI(ns("module2")),
    bsModal(ns("modalhelpvalue"), "Help for p-value distribution", ns("help_pvalue"),size = "large",
            helpText("Possible histogram versions and their meaning"),
            helpText("Case 1:"),
            fluidRow(column(5,img(src="p-value1.PNG"))),
            helpText("If your p-value distribution looks like the above then
                     You have (on the surface) a set of well-behaved p-values.
                     That flat distribution along the bottom is all your null p-values, which are uniformly distributed between 0 and 1.
                     Notice a rectangular shape with a peak at 0 .The taller the peak, the more p-values are close to 0 and therefore significant."),
            helpText("Case 2:"),
            fluidRow(column(5,img(src="p-value2.png"))),
            helpText("The above distribution is bimodal since we have two peaks.This u shaped histogram
                     means low variance and therefore we need to correct for it.To correct for it please close this tab and
                     select the checkboxes corresponding to the p-values to be corrected.
                     "),
            helpText("Case 3:"),
            fluidRow(column(5,img(src="p-value3.png"))),
            helpText("The problem with this distribution is that the N(0,1) null distribution of the Wald
                     test was not appropriate for the dataset. In particular, the assumed variance of the
                     null distribution was too high, and therefore we see a hill-shaped histogram of p-values
                     distribution .To correct for it please close this tab and
                     select the checkboxes corresponding to the p-values to be corrected."),
            helpText("Case 4:"),
            fluidRow(column(5,img(src="p-value4.png"))),
            helpText("The above is a flat distribution (what statisticians call a uniform distribution).
                     At most a small percentage of hypotheses are non-null.No need to correct for such a distribution.")

            ),
    bsModal(ns("modalqc"), "Quality Control", ns("p-value"),size = "large",
            tabPanel("Distribution of p-values plot",
                                      fluidPage(
                                        bscols(widths = 10,div(style="height: 15px;",width = '200px')),
                                        helpText("A histogram of p-value distribution lets you get an immediate sense of how your test
                                                 behaved across all your hypotheses, and helps you to immediately diagnose
                                                 some potential problems.Therefore please have a look at the p-value
                                                 histogram for all selected comparisons. To see how a correct p value distribution
                                                 should look like please click the button below."),
                                        actionButton(ns("help_pvalue"),"Info"),

                                        bscols(widths = 10,div(style="height: 15px;",width = '200px')),
                                        p("To correct for p-values"),
                                        p("Select checkboxes that corresponds to comparisons
                                          that need to be corrected."),
                                        p("Close this window and click start."),
                                        bscols(widths = 10,div(style="height: 15px;",width = '200px')),
                                        bscols(
                                          helpText("For more information regarding p-value distribution"),
                                          a("Click here", href="http://varianceexplained.org/statistics/interpreting-pvalue-histogram/", target="_blank")),
                                        uiOutput(ns("p_comb")),
                                        uiOutput(ns("p_val")),
                                        downloadButton(ns('download_p_value_plot'), 'Download Plot'),
                                        plotOutput(ns("p_value_plot")),
                                        actionButton(ns("help"),"Help"))
                                        ))
  )
}
Module_Differential_Expression<-function(input,output,session,conchoice,dds.fc,
                                         CoCena, normal,batch_choice,batch_corrected,anova_table,conchoice_module1,filechoice)

{

  #Step 1: Define comparisons:
  #We place all the treatment/condition groups in two groups namely 1 and 2
  #Treatment/conditions selected by the user in both groups are combined to form comparisons (treatment 1 vs treatment 2)
  #Among all possible comparisons (for instance 1 vs 2 is same as 2 vs 1,only sign of foldchange will be opposite )
  #the user then selects the comparisons of interest

  #First set of condtion/treatment groups
  output$group1<-
    renderUI({
        pheno<-colData(dds.fc()[[1]])
        if(as.numeric(filechoice())!=3)condition<-pheno[,as.numeric(conchoice)]
        else condition<-pheno[,as.numeric(conchoice_module1())]
        selectInput(session$ns("candidates2"),multiple=TRUE,label = h5("Choose candidate levels for group 1 (e.g. control)") ,
                    choices = condition)
    })

  #Second set of condition/treatment groups
  output$group2<-
    renderUI({
        print("filechoice")
        print(filechoice())
        pheno<-colData(dds.fc()[[1]])
        if(as.numeric(filechoice())!=3)condition<-pheno[,as.numeric(conchoice)]
        else condition<-pheno[,as.numeric(conchoice_module1())]
        selectInput(session$ns("candidates1"),multiple=TRUE,label = h5("Choose candidate levels for group 2 (e.g. treatment)") ,
                    choices = condition)
    })

  #Combination of treatment/condition groups selected in group 1 and group 2
  #Let's suppose group 1 and group 2 have treatment/conditions groups A,B and C
  #This reactive indicates the comparisons user wants for instance A vs B, C vs B .etc
  combination<-reactive({
    req(input$combo)
    if(length(input$candidates1)>0 && length(input$candidates2)>0)
    {
      l<-list()
      for( i in 1:length(input$candidates1))
      {
        for(j in 1:length(input$candidates2))
        {
          if(input$candidates1[i]!=input$candidates2[j])
          {
            if(!(c(input$candidates2[j],input$candidates1[i]) %in% l))
              list2 <- list()
            list2[1] = input$candidates1[i]
            list2[2] = input$candidates2[j]
            l[[length(l)+1]]<-list2
          }
        }
      }

      l
    }
  })
  #Once the above step is computed (all possible combinations of treatment/condition groups between
  #groups 1 and 2), the list is provided as input choices to user.
  #choose the comparisons

                 output$group3<-
                   renderUI({
                     if(input$combo > 0){
                     if(!is.null(combination())){
                       num <- length(combination())
                       choices<-lapply(1:num, function(i) {
                         paste(combination()[[i]][1],' vs ',combination()[[i]][2])

                       })
                       selectInput(session$ns("combination"),multiple=TRUE,label = h5("Choose comparisons") ,
                                                     choices = choices)
                     }
                     }
                   })


  #Once the desired comparisons are selected by the user,we require the following values as input to
  #test for differential expression. The required inputs are
  #1) P-value cut-off
  #2) Fold change cut-off
  #3) Choice of hypithesis (For more information please refer to DESeq2 vignette)
  #4) FDR correction needed? (Following hypothesis testing, we do multiple testing to account for false positive
                              #However in some cases, the number of differentially expressed genes
                              #identified is very low due to multiple testing. Switiching off the
                              #multiple testing can increase the numbe of genes. However these contain false postivies
                              #and is not recommended.)

  #Input 1: define slider for p-value
  output$sliders <- renderUI({
    if(!is.null(input$combination))
    {
      num <- length(input$combination)
      lapply(1:num, function(i) {
        div(style="height: 85px;",sliderInput(inputId = session$ns(paste0("p-value", i)), label = paste("p-value cut-off"),
                                              min = 0.01, max = 0.99, value =0.05, step = 0.01,width = '200px'))

      })
    }
  })

  #Input 2: define sliders for fold change
  output$sliders_fc <- renderUI({

      req(input$combination)
      num <- length(input$combination)
      lapply(1:num, function(i) {
        req(input[[paste0("hypothesis_choice", i)]])
        if(as.numeric(input[[paste0("hypothesis_choice", i)]])==1)
        {
          div(style="height: 85px;",width = '200px')
        }
        else if(as.numeric(input[[paste0("hypothesis_choice", i)]]) %in% c(2,4,5)){
          div(style="height: 85px;",sliderInput(inputId = session$ns(paste0("FC_cutoff", i)), label = paste("FC cut-off"),
                                                min=0.25,max=20,value = 2, step = 0.25,width = '200px'))
        }
        else if(as.numeric(input[[paste0("hypothesis_choice", i)]])==3)
        {
          div(style="height: 85px;",sliderInput(inputId = session$ns(paste0("FC_cutoff", i)), label = paste("FC cut-off"),
                                                min=0.1,max=0.99,value = 0.5, step = 0.1,width = '200px'))
        }

      })
  })

  #Input 3: options for hypothesis choices
  #please refer to DESEq2 manual for more information on hypothesis choices. Use keyword "althernative hypothesis"
  output$choices <- renderUI({

      req(input$combination)
      num <- length(input$combination)
      lapply(1:num, function(i) {
        div(style="height: 95px;", radioButtons(inputId = session$ns(paste0("hypothesis_choice",i)), label = "Choose the hypothesis",
                                                choices = list("default"=1,"greaterAbs"=2,"lessAbs"=3,"greater"=4,"less"=5),
                                                selected = 1,inline = T,width = '400px'))
        })
  })

  #Input 4: Radio button: FDR correction needed?, deafult setting is True(multiple testing is carried out)
  output$ind <- renderUI({

    req(input$combination)
      num <- length(input$combination)
      lapply(1:num, function(i) {
        div(style="height: 95px;", radioButtons(inputId = session$ns(paste0("ind_choice",i)), label = "FDR correction",
                                                choices = list("Yes"=TRUE,"No"=FALSE),
                                                selected = TRUE,inline = T,width = '200px'))

      })
  })
  
  # Display combinations
  output$comb <- renderUI({
      req(input$combination)
      num <- length(input$combination)
      lapply(1:num, function(i) {
        div(style="height: 85px;",p(width='20px',h5(input$combination[i])))

      })
  })

#Now we test for differential expression using the inputs obtained in the previous steps.
#The DESeq function from DESeq2 package is used to test for differential expression. This
#method uses the statistical test "Wald test" to test for differential expression.
#Differential expression analysis=hypothesis testing + Quality control + multiple testing
#Quality control=Check histogram of p-values. If it is not anti-conservative then we correct for it and
#then do multiple testing.
#This reactive automaticaly corrects for p values if p-value hisogram is incorrect

###Input provided
#Comparison eg. A vs B
#P value cut off
#hypothesis
# fold change cut off
######Output returned
#returns a list conatining (up regulated genes, down regulated genes,all genes)

  DE_genes<-eventReactive(input$ok3,{

    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = 'Processing data',
                 value = 0)
    num <- length(input$combination)
    n<-num+2


    library('fdrtool')
    #Get the deseq2 dataset
    dds<-dds.fc()[[1]]
    dds.fc <- dds
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 1,"/",n))

    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)

    #Get the hypothesis selected for each comparison
    hyp_choice<-lapply(1:num, function(i) {
      as.numeric(input[[paste0("hypothesis_choice", i)]])
    })
    #Map hypothesis id to hyp vector
    hyp<-list("default","greaterAbs","lessAbs","greater","less")
    #multicore
    ddsNoPrior<-NULL
    if(3 %in% hyp_choice) ddsNoPrior<-DESeq(dds.fc,betaPrior = FALSE)
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 2,"/",n))

    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    #Get base mean of each condition
    base_mean<-sapply( levels(dds$condition), function(lvl) rowMeans( normal()[,dds$condition == lvl] ) )
    result<-list()

    #additional DESeq
    dds <- DESeq(dds)

    for (i in 1:num)
    {
      cond1<-as.character(strsplit(input$combination[i],"  vs  ",fixed = T)[[1]][1])
       print("cond1")
       print(cond1)
      cond2<-as.character(strsplit(input$combination[i],"  vs  ",fixed = T)[[1]][2])
       print("cond2")
       print(cond2)

      if(hyp_choice[[i]]==1) # if hypothesis is default (null hypothesis is: there is no difference between any two condition/treatment groups)
      {
        res<-NULL
        if(length(levels(dds.fc$condition))==2) #If there are only two treatment/condition groups in an experiment(for eg. control vs knock out)
        {
          if(input[[paste0("ind_choice",i)]]==FALSE){ #If FDR option is selected as false (no multiple testeing is done)
            res <- results(dds,pAdjustMethod = "none", contrast=c("condition",cond1,cond2))
            print("res")
            print(res)
            res.shr <- lfcShrink(dds=dds, contrast=c("condition",cond1,cond2), res=res)
          }


          else  {
            res <- results(dds,pAdjustMethod = "BH", contrast=c("condition",cond1,cond2))
            res.shr <- lfcShrink(dds=dds,type = "normal", contrast=c("condition",cond1,cond2), res=res)

          }
        }
        else{ #If there are more than two treatment/condition groups in an experiment(for eg. control vs knock out)
          if(input[[paste0("ind_choice",i)]]==FALSE){ #If FDR option is selected as false (no multiple testeing is done)
            res <- results(dds,pAdjustMethod = "none",contrast=c("condition",cond1,cond2))
            res.shr <- lfcShrink(dds=dds,type = "normal",contrast=c("condition",cond1,cond2))
          }
          else{
            res <- results(dds,pAdjustMethod = "BH", contrast=c("condition",cond1,cond2))
            res.shr <- lfcShrink(dds=dds,type = "normal", contrast=c("condition",cond1,cond2), res=res)
            }
        }

        if (as.numeric(input$shrinkornot == 2)){
          ma <- res
          res <- res
        }
        else if (as.numeric(input$shrinkornot == 1)){
          ma <- res.shr
          res <- res.shr
        }

        summary(res)
        fc<-lapply(res$log2FoldChange,function(x)
          if(sign(x)==1)  2^x
          else if(sign(x)==-1)((-1) *(1/(1/(2^abs(x)))))
          else if(is.na(x)) NA
        )
        na<-which(is.na(fc))
        non_na<-which(!is.na(fc))
        res$foldchange[non_na]<-as.numeric(unlist(fc[non_na]))
        res$foldchange[na]<-fc[na]
        colnames(res)[7]<-"FoldChange"
        c<-colnames(res)
        r<-rownames
        res <- res[,c(c[1],c[7],c[2],c[3],c[4],c[5],c[6])]
        #p value correction block
        #User selects the comparisons whose P-value need to be corrected. The index is contained in input$p_val
        if(i %in% as.numeric(input$p_val))
        {

          hist(res[,"padj"], col = "royalblue4",
               main = "WT-med vs KO-med initial null model", xlab = "CORRECTED p-values")
          #We first remove genes filtered out by independent filtering and the dispersion outliers,
          #they have NA adj. pvals and NA p-values respectively.
          res <- res[ !is.na(res$padj), ]
          ma <- ma[!is.na(ma$padj),]
          res <- res[ !is.na(res$pvalue), ]

          ma <-ma[ !is.na(ma$pvalue), ]
          #We now remove the original adjusted p-values, since we will add the corrected ones later on.
          res <- res[, -which(names(res) == "padj")]
          ma <- ma[, -which(names(ma) == "padj")]
          FDR.DESeq2Res <- fdrtool(res$stat, statistic= "normal")
          FDR.DESeq2ma <- fdrtool(ma$stat, statistic= "normal")
          FDR.DESeq2Res$param[1, "sd"]
          FDR.DESeq2ma$param[1, "sd"]
          res[,"pvalue"]<-FDR.DESeq2Res$pval
          ma[,"pvalue"]<-FDR.DESeq2ma$pval
          hist(FDR.DESeq2Res$pval, col = "royalblue4",
               main = "WT-med vs KO-med correct p-val null model", xlab = "CORRECTED p-values")
          print("prior to multiple testing")
          #multiple testing on corrected p-values
          if(input[[paste0("ind_choice",i)]]==FALSE) #If user has set FDR to false then multiple testing is waived off
          {
            res[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            ma[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }
          else
          {
            print("inside multiple testing")
            res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
            ma[,"padj"]  <- p.adjust(FDR.DESeq2ma$pval, method = "BH")
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }

        }

        #prepare output
        #Remove genes with foldhcange values =0
        resSig<- subset(res,log2FoldChange <0 || log2FoldChange>0)
        print(nrow(resSig))

        #Merge the basemeans with the DE table
        temp1<-merge(resSig,base_mean,by=0,all=TRUE)
        #filter DE table most significant genes with fdr cut off 0.05
        temp2 <- subset(temp1, pvalue < input[[paste0("p-value",i)]])
        if(input[[paste0("ind_choice",i)]]==TRUE)
          #filter DE table most significant genes with fdr cut off 0.05
          temp2 <- subset(temp1, padj < input[[paste0("p-value",i)]])

        #split temp2 into up and down regulated genes

        #up regulated genes
        up_genes<-subset(temp2,log2FoldChange>0)

        #down regulated genes
        down_genes<-subset(temp2,log2FoldChange<0)

        result[[length(result)+1]]<-list(up_genes,down_genes,temp2,res,ma)
        #output(up regulated genes, down regualted genes, all genes, output table returned by DESeq2 and input for MA plot)

      }
      else if(hyp_choice[[i]]==3) # if hypothesis is lesserAbs(aim is to look for weakly expressed genes)
      {
        res<-NULL
        if(length(levels(dds.fc$condition))==2){
          if(input[[paste0("ind_choice",i)]]==FALSE){
            res <- results(ddsNoPrior,pAdjustMethod = "none",
                           altHypothesis=hyp[hyp_choice[[i]]],
                           lfcThreshold = log2(1/input[[paste0("FC_cutoff",i)]]))
            res.shr <- lfcShrink(dds=ddsNoPrior, contrast=c("condition",cond1,cond2), res=res)}

          else {res <- results(ddsNoPrior,pAdjustMethod = "BH",
                              altHypothesis=hyp[hyp_choice[[i]]],
                              lfcThreshold = log2(1/input[[paste0("FC_cutoff",i)]]))
                res.shr <- lfcShrink(dds=ddsNoPrior, contrast=c("condition",cond1,cond2), res=res)}
        }
        else{
          if(input[[paste0("ind_choice",i)]]==FALSE) {#if FDR is set to false by user then multiple testing is waived off
            res <- results(ddsNoPrior,pAdjustMethod = "none",
                           altHypothesis="lessAbs",
                           contrast=c("condition",cond1,cond2),
                           lfcThreshold = log2(1/input[[paste0("FC_cutoff",i)]]))
            res.shr <- lfcShrink(dds=ddsNoPrior, contrast=c("condition",cond1,cond2), res=res)}
          else{ res <- results(ddsNoPrior,pAdjustMethod = "BH",
                              altHypothesis="lessAbs",
                              contrast=c("condition",cond1,cond2),
                              lfcThreshold = log2(1/input[[paste0("FC_cutoff",i)]]))
                res.shr <- lfcShrink(dds=ddsNoPrior, contrast=c("condition",cond1,cond2), res=res)}
        }

        if (as.numeric(input$shrinkornot == 2)){
          ma <- res
          res <- res
        }
        else if (as.numeric(input$shrinkornot == 1)){
          ma <- res.shr
          res <- res.shr
        }

        fc<-lapply(res$log2FoldChange,function(x)
          if(sign(x)==1)  2^x
          else if(sign(x)==-1) 1/((-1) *(1/(1/(2^abs(x)))))
          else if(is.na(x)) NA
        )

        na<-which(is.na(fc))
        non_na<-which(!is.na(fc))
        res$foldchange[non_na]<-as.numeric(fc[non_na])
        res$foldchange[na]<-fc[na]
        colnames(res)[7]<-"FoldChange"
        c<-colnames(res)
        res <- res[,c(c[1],c[7],c[2],c[3],c[4],c[5],c[6])]

        #Quality control. User selectes the comparisons whose p-value needs to be corrected.
        #input$p_val contains the index of the selected comparisons
        if(i %in% as.numeric(input$p_val))
        {
          #We first remove genes filtered out by independent filtering and the dispersion outliers,
          #they have NA adj. pvals and NA p-values respectively.
          res <- res[ !is.na(res$padj), ]
          ma <- ma[ !is.na(ma$padj), ]
          res <- res[ !is.na(res$pvalue), ]
          ma <- ma[ !is.na(ma$pvalue), ]

          #We now remove the original adjusted p-values, since we will add the corrected ones later on.
          res <- res[, -which(names(res) == "padj")]
          ma <- ma[, -which(names(ma) == "padj")]
          FDR.DESeq2Res <- fdrtool(res$stat, statistic= "normal")
          FDR.DESeq2ma <- fdrtool(ma$stat, statistic= "normal")
          FDR.DESeq2Res$param[1, "sd"]
          FDR.DESeq2ma$param[1, "sd"]
          res[,"pvalue"]<-FDR.DESeq2Res$pval
          ma[,"pvalue"]<-FDR.DESeq2ma$pval

          #multiple testing on corrected p-values
          if(input[[paste0("ind_choice",i)]]==FALSE) #If user has set FDR to false then multiple testing is waived off
          {
            res[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            ma[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }
          else
          {
            res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
            ma[,"padj"]  <- p.adjust(FDR.DESeq2ma$pval, method = "BH")
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }
        }
        #Merge the basemeans with the DE table
        temp<-merge(res,base_mean,by=0,all=TRUE)
        #filter DE table most significant genes with fdr cut off 0.05
        resSig <- subset(temp, pvalue < input[[paste0("p-value",i)]])
        if(input[[paste0("ind_choice",i)]]==TRUE)
          resSig <- subset(temp, padj < input[[paste0("p-value",i)]])
        #up regulated genes
        up_genes<-subset(resSig,log2FoldChange>0)
        #down regulated genes
        down_genes<-subset(resSig,log2FoldChange<0)

        result[[length(result)+1]]<-list(up_genes,down_genes,resSig,res,ma)
        # return(up regulated genes,down regulated genes, all genes, output returned by DESeq, input for ma plot)

      }

      else if(hyp_choice[[i]] %in% c(2,4,5)){ #if the hypothesis is greaterAbs(consider only genes with logfoldchange > threshold as significantly differentially expressed)
                                                                    #greater(aim to look for up regulated genes), if maximum of expected differentially expressed genes is up regulated
        res<-NULL                                                   #lesser(aim is to look for down regulated genes), if maximum of expected differentially expressed genes is downregulated
        if(length(levels(dds.fc$condition))==2) #If there are only two treatment/condition groups in an experiment
        {
          if(input[[paste0("ind_choice",i)]]==FALSE){
            res <- results(dds,pAdjustMethod = "none", #If FDR is set by user to false, then multiple testing is waived off
                           altHypothesis=hyp[[hyp_choice[[i]]]],
                           lfcThreshold = log2(input[[paste0("FC_cutoff",i)]]))#,
            res.shr <- lfcShrink(dds=dds, contrast=c("condition",cond1,cond2), res=res)}

          else{ res <- results(dds,pAdjustMethod = "BH",
                              altHypothesis=hyp[[hyp_choice[[i]]]],
                              lfcThreshold = log2(input[[paste0("FC_cutoff",i)]]))
          res.shr <- lfcShrink(dds=dds, contrast=c("condition",cond1,cond2), res=res)}

        }
        else{#If there are more than two treatment/condition groups in an experiment
          if(input[[paste0("ind_choice",i)]]==FALSE){
            res <- results(dds,pAdjustMethod = "none",#If FDR is set by user to false, then multiple testing is waived off
                           altHypothesis=hyp[[hyp_choice[[i]]]],
                           contrast = c("condition",cond1,cond2),
                           lfcThreshold = log2(input[[paste0("FC_cutoff",i)]]))
            res.shr <- lfcShrink(dds=dds, contrast=c("condition",cond1,cond2), res=res)}

          else{ res <- results(dds,pAdjustMethod = "BH",
                              altHypothesis=hyp[[hyp_choice[[i]]]],
                              contrast = c("condition",cond1,cond2),
                              lfcThreshold = log2(input[[paste0("FC_cutoff",i)]]))
                res.shr <- lfcShrink(dds=dds, contrast=c("condition",cond1,cond2), res=res)}

        }

        if (as.numeric(input$shrinkornot == 2)){
          ma <- res
          res <- res
        }
        else if (as.numeric(input$shrinkornot == 1)){
          ma <- res.shr
          res <- res.shr
        }

        fc<-lapply(res$log2FoldChange,function(x)
          if(sign(x)==1)  2^x
          else if(sign(x)==-1) 1/((-1) *(1/(1/(2^abs(x)))))
          else if(is.na(x)) NA
        )

        na<-which(is.na(fc))
        non_na<-which(!is.na(fc))
        res$foldchange[non_na]<-as.numeric(fc[non_na])
        res$foldchange[na]<-fc[na]
        colnames(res)[7]<-"FoldChange"
        c<-colnames(res)
        res <- res[,c(c[1],c[7],c[2],c[3],c[4],c[5],c[6])]
        if(i %in% as.numeric(input$p_val))
        {
          #We first remove genes filtered out by independent filtering and the dispersion outliers,
          #they have NA adj. pvals and NA p-values respectively.
          res <- res[ !is.na(res$padj), ]
          ma <- ma[ !is.na(ma$padj), ]
          res <- res[ !is.na(res$pvalue), ]
          ma <- ma[ !is.na(ma$pvalue), ]
          #We now remove the original adjusted p-values, since we will add the corrected ones later on.
          res <- res[, -which(names(res) == "padj")]
          ma <- ma[, -which(names(ma) == "padj")]
          FDR.DESeq2Res <- fdrtool(res$stat, statistic= "normal")
          FDR.DESeq2ma <- fdrtool(ma$stat, statistic= "normal")
          FDR.DESeq2Res$param[1, "sd"]
          FDR.DESeq2ma$param[1, "sd"]
          res[,"pvalue"]<-FDR.DESeq2Res$pval
          ma[,"pvalue"]<-FDR.DESeq2ma$pval

          #multiple testing on corrected p-values
          if(input[[paste0("ind_choice",i)]]==FALSE) #If user has set FDR to false then multiple testing is waived off
          {
            res[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            ma[,"padj"]  <- FDR.DESeq2Res$pval #p.adj has the same value as p value
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }
          else
          {
            res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
            ma[,"padj"]  <- p.adjust(FDR.DESeq2ma$pval, method = "BH")
            hist(res[,"padj"], col = "royalblue4",
                 main = "WT-med vs KO-med correctp-adj  null model", xlab = "CORRECTED p-values")
          }
        }
        #Merge the basemeans with the DE table
        temp<-merge(res,base_mean,by=0,all=TRUE)
        #filter DE table most significant genes with fdr cut off 0.05
        resSig <- subset(temp, pvalue < input[[paste0("p-value",i)]])
        if(input[[paste0("ind_choice",i)]]==TRUE)
          resSig <- subset(temp, padj < input[[paste0("p-value",i)]])
        #up regulated genes
        up_genes<-subset(resSig,log2FoldChange>0)
        #down regulated genes
        down_genes<-subset(resSig,log2FoldChange<0)

        result[[length(result)+1]]<-list(up_genes,down_genes,resSig,res,ma)

      }


      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", i+2,"/",n))

      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }

    if(!is.null(CoCena()))
    {
      
      cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
        # filter anova table from module genes
      for(i in 1:nrow(cluster_info)){
            gene_cluster<-unlist(strsplit(cluster_info[i,"gene_n"],split=","))
          result[[length(result)+1]]<-list(0,0,gene_cluster,0,0)
       
          # Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = paste("Doing part", num+2+i,"/",n))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)

      }
    }

    result
  })

  #display p-value plot
  observeEvent(input$ok3, {
    toggleModal(session, "modalqc", toggle = "open")
  })

  #Display summary of DE gene table
  output$de_genes <- DT::renderDataTable({
      req(input$ok3)
    if(input$ok3>0)
    {
      print("line 728")
      result<-DE_genes()

      num <- length(input$combination)
      rows<-num
      res<-data.frame(matrix(NA, nrow = num, ncol = 3))

      if(!is.null(CoCena()))
      {
        cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
        print("cluster_info")
        print(cluster_info)
        res<-data.frame(matrix(NA, nrow = num+nrow(cluster_info), ncol = 3))
        colnames(res)<-c('Up regulated','Down regulated','Both')
        rowname_vector<-c()
          for(i in 1:num)
          {
            res[i,1]<-nrow(as.data.frame(result[[i]][1]))
            res[i,2]<-nrow(as.data.frame(result[[i]][2]))
            res[i,3]<-nrow(as.data.frame(result[[i]][3]))
            rowname_vector[i]<-input$combination[i]
          }
          for(k in (num+1):(nrow(cluster_info)+num))
          {
            res[k,1:2]<-0
            res[k,3]<-length(result[[k]][[3]])
            rowname_vector[k]<-cluster_info$color[k-num]
          }
          rownames(res)<-rowname_vector
        
      }

      else{
        rownames(res)<-lapply(1:num, function(i) {
          input$combination[i]
        })
        colnames(res)<-c('Up regulated','Down regulated','Both')
        for(i in 1:num)
        {
          res[i,1]<-nrow(as.data.frame(result[[i]][1]))
          res[i,2]<-nrow(as.data.frame(result[[i]][2]))
          res[i,3]<-nrow(as.data.frame(result[[i]][3]))
        }
      }

      DT::datatable(res,class = 'cell-border stripe',
                    selection = list(mode='single',target = 'cell'),
                    extensions = list('Scroller'=NULL,'Buttons'=NULL),
                    options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                   buttons = list(list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'), text = 'Download table'))),#I('colvis')

                    escape = FALSE)

    }

  })



  #on clicking a specific comparison, the corresponding table(list of DE genes) should appear
  #On clicking a module the anova table is displayed
  observeEvent(input$de_genes_cell_clicked,{
    selected <- input$de_genes_cells_selected
    row<-selected[1]
    col<-selected[2]
    if(length(selected)>0){

      #table diaplyas the DE genes
      output$filtered_data <- DT::renderDataTable({
        if(!is.null(CoCena()))
        {
          cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
            num <- length(input$combination)
            if(((col==3) && (row>length(input$combination))))
            {
            
              gene_cluster<-unlist(strsplit(cluster_info[(row-num),"gene_n"],split=","))

              #########preparing the anova table in the order as output#########
              a_tab<-anova_table()[,-c(2,3)]
              cond<-unique(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
              c<-colnames(a_tab)
              temp<-as.vector(c[4:(3+length(cond))])
              temp2<-as.vector(c[(length(cond)+4):length(c)])
              all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
              ##################################################################
              anova_genes<-rownames(all_genes)
              df<-all_genes[which(anova_genes %in% gene_cluster),]

              DT::datatable(df,class = 'cell-border stripe',
                            selection = list(mode='single',target = 'row'),
                            extensions = list('Scroller'=NULL,'Buttons'=NULL),
                            options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                           buttons = list()))

            }
            else if(((col<3) && (row>length(input$combination)))){
              df<-data.frame(gene_list=character())
              DT::datatable(df,class = 'cell-border stripe',
                            selection = list(mode='single',target = 'row'),
                            extensions = list('Scroller'=NULL,'Buttons'=NULL),
                            options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                           buttons = list()))
            }
            else if(row<=length(input$combination))
            {
              result<-DE_genes()
              df<-as.data.frame(result[[row]][col])
              genes<-df[,1]
              df<-df[-1]
              rownames(df)<-genes
              df_final<- df[order(df$padj),]
              colnames(df_final)[1]<-"Overall mean"
              for (i in 8:length(colnames(df_final)))
              {
                colnames(df_final)[i]<-paste(colnames(df_final)[i],"mean")
              }
              a_tab<-as.data.frame(df_final[,-4])
              c<-colnames(a_tab)
              temp<-as.vector(c[7:length(c)])
              tab <- a_tab[,c(temp,c[2],c[3],c[5],c[6],c[1],c[4])]
              df_de=NULL
              if(input$lfc_de =="No")
              {
                names<-colnames(tab)
                idx<-which(names =="log2FoldChange")
                df_de=tab[,-idx]
              }
              else df_de=tab
              DT::datatable(df_de,class = 'cell-border stripe',
                            selection = list(mode='single',target = 'row'),
                            extensions = list('Scroller'=NULL,'Buttons'=NULL),
                            options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                           buttons = list()))

            }
          
        }

        else
        {
          result<-DE_genes()
          df<-as.data.frame(result[[row]][col])
          genes<-df[,1]
          df<-df[-1]
          rownames(df)<-genes
          df_final<- df[order(df$padj),]
          colnames(df_final)[1]<-"Overall mean"
          for (i in 8:length(colnames(df_final)))
          {
            colnames(df_final)[i]<-paste(colnames(df_final)[i],"mean")
          }
          a_tab<-as.data.frame(df_final[,-4])
          c<-colnames(a_tab)
          temp<-as.vector(c[7:length(c)])
          tab <- a_tab[,c(temp,c[2],c[3],c[5],c[6],c[1],c[4])]
          df_de=NULL
          if(input$lfc_de =="No")
          {
            names<-colnames(tab)
            idx<-which(names =="log2FoldChange")
            df_de=tab[,-idx]
          }
          else df_de=tab
          DT::datatable(df_de,class = 'cell-border stripe',
                        selection = list(mode='single',target = 'row'),
                        extensions = list('Scroller'=NULL,'Buttons'=NULL),
                        options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                       buttons = list()))

        }

      })

      ###plotcounts of gene
      #display gene expression across conditon
      observeEvent(input$filtered_data_rows_selected,{
        #get which gene was clicked
        #get the row  number of the gene clicked
        selected_row <- input$filtered_data_rows_selected
        #get the normalized data
        full_data<-normal()
        #get the DE table where we clicked the gene
        selected <- input$de_genes_cells_selected
        row<-selected[1]
        col<-selected[2]
        if(length(selected)>0)
        {
          #get the table where we clicked the gene
          result<-DE_genes()
          df<-NULL
          df_final<-NULL
          an_gene<-NULL
          if(row<=length(input$combination))
          {
            df<-as.data.frame(result[[row]][col])
            df_final<- df[order(df$padj),]
            print(head(df_final))
            an_gene<-df_final[selected_row,1]
          }
          else if(row>length(input$combination))
              {
                cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
                num <- length(input$combination)
                gene_cluster<-unlist(strsplit(cluster_info[(row-num),"gene_n"],split=","))
                
                #########preparing the anova table in the order as output#########
                a_tab<-anova_table()[,-c(2,3)]
                cond<-unique(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
                c<-colnames(a_tab)
                temp<-as.vector(c[4:(3+length(cond))])
                temp2<-as.vector(c[(length(cond)+4):length(c)])
                all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
                ##################################################################
                print(head(all_genes))
                anova_genes<-rownames(all_genes)
                df_final<-all_genes[which(anova_genes %in% gene_cluster),]
                an_gene<-rownames(df_final)[selected_row]
              }

          counts<-as.vector(full_data[which(an_gene %in% rownames(full_data)),])
          if (as.numeric(filechoice())==3) condition <-cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice_module1())]) #already chosen for STAR
          else cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
          library('data.table')
          df<-data.frame(counts,cond)
          colnames(df)<-c('count','condition')
          batch<-batch_choice()
          if(as.numeric(batch)==1) callModule(gene_count_module,"module2",NULL,reactive({dds.fc()[[1]]}),reactive({an_gene}))
          else
            {
              idx<-which(rownames(batch_corrected()) %in% an_gene)
              callModule(gene_count_module,"module2",batch_corrected()[idx,],reactive({dds.fc()[[1]]}),reactive({an_gene}))
            }
        }

        })
      #download button
      output$download_DE_genes_Table <- downloadHandler(

        filename = function()
        {
          if(row<=length(input$combination))
          {
            condition<-input$combination[row]
            if(col==1 & as.numeric(input$datachoice4==1)) paste('Up regulated genes for ',condition,'.xlsx')
            else if (col==1 & as.numeric(input$datachoice4==2)) paste('Up regulated genes for ',condition,'.csv')
            else if(col==2 & as.numeric(input$datachoice4==1)) paste('Down regulated genes for ',condition,'.xlsx')
            else if(col==2 & as.numeric(input$datachoice4==2)) paste('Down regulated genes for ',condition,'.csv')
            else if (as.numeric(input$datachoice4==1)) paste('DE genes for ',condition,'.xlsx')
            else if (as.numeric(input$datachoice4==2)) paste('DE genes for ',condition,'.csv')
          }
          else{
            if (as.numeric(input$datachoice4==1)){
              cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
              num <- length(input$combination)
              paste('ANOVA genes for ',cluster_info$color[row-num],'.xlsx')
            }
            else if (as.numeric(input$datachoice4==2)){
              cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
              num <- length(input$combination)
              paste('ANOVA genes for ',cluster_info$color[row-num],'.csv')
            }

          }

        },
        content = function(file) {
        if (as.numeric(input$datachoice4==1)){
          nam<-"Sheet 1"
          df_de<-NULL
          if(row<=length(input$combination))
          {
            #sort by adjusted p value.
            result<-DE_genes()
            df<-as.data.frame(result[[row]][col])
            genes<-df[,1]
            df<-df[-1]
            rownames(df)<-genes
            df_final<- df[order(df$padj),]
            colnames(df_final)[1]<-"Overall mean"
            for (i in 8:length(colnames(df_final)))
            {
              colnames(df_final)[i]<-paste(colnames(df_final)[i],"mean")
            }
            a_tab<-as.data.frame(df_final[,-4])
            print(colnames(a_tab))
            c<-colnames(a_tab)
            print(length(c))
            temp<-as.vector(c[7:length(c)])
            df_de <- a_tab[,c(temp,c[2],c[3],c[5],c[6],c[1],c[4])]
            condition<-input$combination[row]
            condition<-str_replace_all(condition,"[^[:alnum:]]",".")
            if(col==1) nam<-paste('Up regulated genes for ',condition)
            else if(col==2) nam<-paste('Down regulated genes for ',condition)
            else nam<-paste('DE genes for ',condition)
          }
          else
          {
            if(row>length(input$combination))
            {
              cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
              num <- length(input$combination)
              gene_cluster<-unlist(strsplit(cluster_info[(row-num),"gene_n"],split=","))
              nam<-paste('ANOVA genes for ',modules$Var1[row-num])
             
              #########preparing the anova table in the order as output#########
              a_tab<-anova_table()[,-c(2,3)]
              cond<-unique(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
              c<-colnames(a_tab)
              temp<-as.vector(c[4:(3+length(cond))])
              temp2<-as.vector(c[(length(cond)+4):length(c)])
              print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
              all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
              ##################################################################
              anova_genes<-rownames(all_genes)
              df_de<-all_genes[which(anova_genes %in% gene_cluster),]
            }
          }

          M <- as.matrix(df_de)
          wb <- createWorkbook()
          addWorksheet(wb, sheetName = "DE Genes")
          writeData(wb = wb, sheet = 1, x = M, colNames = T, rowNames = T)
          saveWorkbook(wb, file)

          }
          else{
            nam<-"Sheet 1"
            df_de<-NULL
            if(row<=length(input$combination))
            {
              #sort by adjusted p value.
              result<-DE_genes()
              df<-as.data.frame(result[[row]][col])
              genes<-df[,1]
              df<-df[-1]
              rownames(df)<-genes
              df_final<- df[order(df$padj),]
              colnames(df_final)[1]<-"Overall mean"
              for (i in 8:length(colnames(df_final)))
              {
                colnames(df_final)[i]<-paste(colnames(df_final)[i],"mean")
              }
              a_tab<-as.data.frame(df_final[,-4])
              c<-colnames(a_tab)
              temp<-as.vector(c[7:length(c)])
              df_de <- a_tab[,c(temp,c[2],c[3],c[5],c[6],c[1],c[4])]

              condition<-input$combination[row]
              condition<-str_replace_all(condition,"[^[:alnum:]]",".")
              if(col==1) nam<-paste('Up regulated genes for ',condition)
              else if(col==2) nam<-paste('Down regulated genes for ',condition)
              else nam<-paste('DE genes for ',condition)
            }
            else
            {
              if(row>length(input$combination))
              {

                cluster_info<-as.data.frame(CoCena()$cluster_calc()[CoCena()$cluster_calc()$cluster_included=="yes",])
                num <- length(input$combination)
                gene_cluster<-unlist(strsplit(cluster_info[(row-num),"gene_n"],split=","))
                nam<-paste('ANOVA genes for ',modules$Var1[row])
             
                #########preparing the anova table in the order as output#########
                a_tab<-anova_table()[,-c(2,3)]
                cond<-unique(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
                c<-colnames(a_tab)
                temp<-as.vector(c[4:(3+length(cond))])
                temp2<-as.vector(c[(length(cond)+4):length(c)])
                all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
                ##################################################################
                anova_genes<-rownames(all_genes)
                df_de<-all_genes[which(anova_genes %in% gene_list),]
              }
            }
            write.csv(df_de, file)}
        }

      )


    }
  })

#######################p-value plot
      #Create drop down menu containing combinations for p value plot( A vs B, C vs D etc) for user selection
      #Selection helps one view p-value histogram for a specific combination
      output$p_comb <- renderUI({
        if(!is.null(combination()))
        {

          num<- length(input$combination)
          comb<-lapply(1:num, function(i) {
            input$combination[i]
          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          selectInput(session$ns("p_choice"),label = h5("Choose comparison") ,
                      choices = checklist,selected = 1)
        }
      })

      #Display check box of all combinations.
      #User chooses combinations for which the p-value needs to be corrected after
      #visualizing p-value histogram of each combination

      output$p_val <- renderUI({
        if(!is.null(combination() ))
        {

          num<- length(input$combination)
          comb<-lapply(1:num, function(i) {
            input$combination[i]

          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          checkboxGroupInput(session$ns("p_val"), "Choose the p-value to be corrected", checklist,selected = NULL)
        }
      })
      #Update the above checkbox if more combinations are added/ if either or both group1 and group 2 input changes
      observe({
        req(input$ok3)
        input$combination
        input$combo
        input$p_val

        if(input$ok3>0 && !is.null(input$combination))
        {
          num<- length(input$combination)
          comb<-lapply(1:num, function(i) {
            input$combination[i]

          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          updateCheckboxGroupInput(session=session, inputId=session$ns("p_val"), choices=checklist, selected=input$p_val)
        }

      })
      observeEvent(input$ok3,
                   {
                     #distribution of p-values
                     output$p_value_plot<- renderPlot({
                       n<-as.numeric(input$p_val)
                       if(!is.null(input$p_choice))
                       {
                         num<- length(input$combination)
                         #Get the deseq2 dataset
                         result<-DE_genes()
                         res<-as.data.frame(result[[as.numeric(input$p_choice)]][4])
                         i<-as.numeric(input$p_choice)
                         hist(res$pvalue[res$baseMean > 1],main = paste(input$combination[i]),  col = "lavender", xlab = "p-values")
                       }


                     })
                     #download p-value plot
                     output$download_p_value_plot <- downloadHandler(
                       filename = paste(" p-value plot.pdf"),
                       content = function(file) {
                         pdf(file)
                         num<- length(input$combination)
                         #Get the deseq2 dataset
                         result<-DE_genes()
                         res<-as.data.frame(result[[as.numeric(input$p_choice)]][4])
                         i<-as.numeric(input$p_choice)

                         hist(res$pvalue[res$baseMean > 1],main = paste("p-value plot for ",input$combination[i]),  col = "lavender", xlab = "p-values")

                         dev.off()
                       })
                   })


########################

return(list(combination=reactive({lapply( 1:length(input$combination), function(i) input$combination[i])}),
            de_genes=reactive({DE_genes()}),
            p_values=reactive({lapply( 1:length(input$combination), function(i) input[[paste0("p-value",i)]] )
            }),
            hypothesis_choice= reactive({lapply( 1:length(input$combination), function(i) input[[paste0("hypothesis_choice",i)]] )
            }),
            ok3=reactive({input$ok3})))


}
