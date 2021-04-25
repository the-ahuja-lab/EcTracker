
### Package installation through CRAN###


install.packages(c("BiocManager","shiny","Seurat","ggplot2","shinythemes","shinyalert","shinyFiles","shinyjs","shinybusy","shinyWidgets","shinyalert","shinyBS","igraph","shinyalert","textshape","tidyr","tidyverse","reshape2","data.table","networkD3","rtracklayer","dplyr","cli","memisc","devtools","memoise"))

### Package installation through Bioconductor###

BiocManager::install(c("TissueEnrich","SummarizedExperiment","SingleCellExperiment","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","karyoploteR","dorothea","bcellViper","viper","AUCell","org.Mm.eg.db","scater","scran","limma","edgeR","dynamicTreeCut","cluster","batchelor","cowplot","DESeq2","MAST", dependencies = T))

####installation through github (Require Rtools)###
library(devtools)
install_github("immunogenomics/harmony")
######## To create EcTracker application###

install.packages("shinyShortcut")
library(shinyShortcut)
shinyShortcut(shinyDirectory = getwd(), OS = .Platform$OS.type,gitIgnore = FALSE)

