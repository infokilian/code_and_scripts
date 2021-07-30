#devtools packages
library(BiocManager)


list.of.dev.packages <- c("plotly",
                          "crosstalk",
                          "DT")


new.packages.dev <- list.of.dev.packages[!(list.of.dev.packages %in% installed.packages()[,"Package"])]

if( "plotly" %in% new.packages.dev) devtools::install_github("ropensci/plotly", upgrade = "never")
if( "crosstalk" %in% new.packages.dev) devtools::install_github("rstudio/crosstalk",force=TRUE, upgrade = "never")
if( "DT" %in% new.packages.dev) devtools::install_github('rstudio/DT', upgrade = "never")
#
#

#CRAN packages
list.of.packages <- c("shiny",
                      "shinyBS",
                      #"shinymaterial",
                      "RcppArmadillo",
                      "shinyjs",
                      'RColorBrewer',
                      "stringr",
                      'formula.tools',
                      'data.table',
                      'fdrtool',
                      "VennDiagram",
                      "devtools",
                      'colorspace',
                      "officer",
                      "magrittr",
                      "openxlsx",
                      "ggrepel",
                      "V8",
                      "WGCNA",
                      'svglite',
                      "visNetwork",
                      "ggpubr",
                      "gplots",
                      "bindrcpp",
                      "pheatmap",
                      "purrr",
                      "flextable",
                      "httr",
                      "jsonlite",
                      "igraph",
                      "MCDA",
                      "ggsci",
                      "gtools",
                      "stringi",
                      "ggplot2"
)

lapply(list.of.packages, require, character.only = TRUE)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages, dependencies = T)

#BioconductoR packages
list.of.bioc.packages<- c("rhdf5",
                          "DESeq2",
                          "IHW",
                          "tximport",
                          'clusterProfiler',
                          "org.Hs.eg.db",
                          "org.Mm.eg.db",
                          "org.Mmu.eg.db",
                          "org.Dm.eg.db",
                          "org.Ce.eg.db",
                          "org.Bt.eg.db",
                          "org.Cf.eg.db",
                          "org.Dr.eg.db",
                          "org.Gg.eg.db",
                          "org.Pt.eg.db",
                          "org.Rn.eg.db",
                          "org.Ss.eg.db",
                          "org.Pf.plasmo.db",
                          "org.At.tair.db",
                          "org.Ag.eg.db",
                          "org.EcK12.eg.db",
                          "org.EcSakai.eg.db",
                          "org.Xl.eg.db",
                          "sva",
                          "limma",
                          "geneplotter",
                          'biomaRt',
                          "AnnotationDbi",
                          "impute",
                          "GO.db",
                          "preprocessCore",
                          "pathview",
                          "lpsymphony",
                          "msigdbr",
                          "ComplexHeatmap"
)


new.packages.bioc <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]

#install.packages("BiocManager")
if(length(new.packages.bioc)>0) BiocManager::install(new.packages.bioc,suppressUpdates=TRUE)

lapply(c(list.of.dev.packages,list.of.packages,list.of.bioc.packages), require, character.only = TRUE)
