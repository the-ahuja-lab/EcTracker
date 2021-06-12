# EcTracker
EcTracker: Tracking and elucidating ectopic expression leveraging large scale scRNA-seq studies

## Webserver 
https://the-ahuja-lab.shinyapps.io/ectracker/
## Workflow
<img src="main/www/workflow.png"> <br/>
###
## Key Points<br/> 
1. EcTracker possesses two distinct modules i.e. CellEnrich and TissueEnrich, to allow detection of cell-type or tissue-specific genes in the user-supplemented scRNA-sequencing dataset.
2. CellEnrich utilizes the widely used AUCell package and leverages on the bonafide fetal and adult tissue markers curated from the large-scale single-cell atlases.
3. EcTracker provides an enhanced visualization method for tracking the chromosomal locations of the selected genes in a cluster-wise manner. 
4. EcTracker implements DoRothEA and its associated statistical method VIPER that allow users to identify central transcription factors driving the expression of ectopic transcripts present in the user datasets.

## Installation of the libraries <br/>
All the libraries can be installed by two main commands: <br/>
1. install.packages(“Package_name”) <br/>
2. install.packages("BiocManager") <br/>
   BiocManager::install("Package_name")<br/>
   
For example:<br/>
To install package karyoploteR, one of the following commands can be used:<br/>
```
install.packages(“karyoploteR”)
# OR
install.packages("BiocManager") 
BiocManager::install("karyoploteR")
```
## List of libraries<br/>
karyoploteR (version: 1.14.1), org.Hs.eg.db (version:3.11.4), TxDb.Hsapiens.UCSC.hg19.knownGene (version:3.2.2), GenomicFeatures (version:1.40.1), GSEABase (version:1.50.1), graph (version:1.66.0), annotate (version:1.66.0), AnnotationDbi (version:1.50.3), GenomicRanges (version:1.40.0), GenomeInfoDb (version:1.24.2), viper (version: 1.22.0), dplyr (version:1.0.3), bcellViper (version:1.24.0), networkD3 (version:0.4), dorothea (version:1.0.1), AUCell (version:1.10.0), shiny (version:1.6.0), shinyFiles (version:0.9.0), shinyalert (version:2.0.0), shinyWidgets (version:0.5.6), shinybusy (version:0.2.2), igraph (version:1.2.6), shinyjs (version:2.0.0), shinycustomloader (version:0.9.0), shinythemes (version:1.2.0), Seurat (version:3.9.9.9008), ggplot2 (version:3.3.3), shinyShortcut <br/>                                                     

## One time quick installation of EcTracker (standalone)<br/>
1. Download the github repository. <br/>
2. User should have latest version of R installed on their system.
3. In R, set the working directory by using the command: <br/>
   ```
   setwd("Ectracker-main/main") 
   ```
4. User should also install Rtools from website https://cran.r-project.org/bin/windows/Rtools/ <br/>
5. Install the dependencies by running line by line code of <b>Dependencies.R</b> on R or RStudio 
6. Once all the dependencies are installed,
```
shinyShortcut(shinyDirectory = getwd(), OS = .Platform$OS.type,gitIgnore = FALSE)
```
This will create a file <b>shinyShortcut.cmd</b> in folder Ectracker-main/main/.shinyrun which the user can use to run the application.<br/>
 
## Instructions to run EcTracker code locally<br/> 
1. Download the github repository. <br/>
2. Install all the required packages. <br/>
3. Put all the files and folders in one directory.<br/>
4. Set the working directory by using the command: <br/>
   ```
   setwd("Ectracker-main/main")
   ```
5. Open the ui.R and server.R in RStudio <br/>
6. To run the app, user can directly click on <b>Run App</b> or use command: 
   ```
   shiny::runApp()
   ```



