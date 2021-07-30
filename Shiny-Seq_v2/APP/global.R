#devtools packages

list.of.dev.packages <- c("plotly",
                          "crosstalk",
                          "DT")

library(plotly)
library(crosstalk)
library(DT)


# CRAN packages
list.of.packages <- c("shiny",
                      "shinyBS",
                      "RcppArmadillo",
                      "GSEABase",
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
                      "purrr"
)

# BioconductoR packages
list.of.bioc.packages<- c("rhdf5",
                          "DESeq2",
                          "IHW",
                          "tximport",
                          'clusterProfiler',
                          "org.Hs.eg.db",
                          "org.Mm.eg.db",
                          "org.Mmu.eg.db",
                          "sva",
                          "limma",
                          "geneplotter",
                          "rhdf5",
                          'biomaRt',
                          "AnnotationDbi", 
                          "impute", 
                          "GO.db", 
                          "preprocessCore",
                          "pcaGoPromoter",
                          "pcaGoPromoter.Mm.mm9",
                          "pcaGoPromoter.Hs.hg19",
                          "pathview",
                          "lpsymphony",
                          "officer"
)

lapply(c(list.of.dev.packages,list.of.packages,list.of.bioc.packages), require, character.only = TRUE)
