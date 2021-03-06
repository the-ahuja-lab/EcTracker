library(shiny)
library(Seurat)
library(shinythemes)
library(shinycustomloader)
library(shinyjs)
library(igraph)
library(shinybusy)
library(XVector)
library(shinyWidgets)
library(shinyalert)
library(shinyFiles)
library(memisc)
library(AUCell)
library(dorothea)
library(networkD3)

library(bcellViper)
library(dplyr)
library(viper)
#library(V8)
jscode <- "shinyjs.toTop = function() {window.scrollTo(0, 0);}" ### for to top button


shinyUI(fluidPage(
  theme=shinytheme("cyborg"),
  tags$style(HTML("
        .tabs-above > .nav > li[class=active] > a {
           background-color: #000;
           color: #FFFF;
        }")),
  
 
 useShinyalert(),
  
  navbarPage(em("EcTracker"),
                   tabPanel("Home" , 
                            icon = icon("home"),
                          
                           HTML('<center><img src="aeae.gif", width="58%", height= "58%" ></center>'),
                         
                           
                          
                          column(8, align="center", offset = 2 ,print("Developed by: "), tags$u("Vishakha Gautam"), print("&" ), tags$u("Siddhant Kalra"),
                          
                          br(),print("For more details:"),
                          tags$head(tags$style(HTML("a {color: blue}"))),
                          tags$a(href="https://ahuja-lab.in/", "The Ahuja Lab", style = "color:yellow"))
                            

                            ),
                   
             
             tabPanel("Analysis" , icon = icon("angle-double-right"),
                     
                     navlistPanel( 
                     tabPanel(tags$b("1. Upload") , 
                               
                                  mainPanel(fileInput("file","Upload Expression Matrix"),
                                            
                                            h5("OR"),br(),
                                            fileInput("File10","Upload 10X Data Files (barcode.tsv, genes.tsv & matrix.mtx) ", multiple = TRUE),
                                            
                                            checkboxInput("testme", " Transdifferentiation (Sample Data)", 
                                                          value = FALSE),

                                            br(),br(),
                                             h6("Information about Data"),
                                           withLoader(tableOutput("dat"),type="image",loader = "ectracker_popup.gif"),
                                           h6("Dimensions"),
                                           withLoader(verbatimTextOutput("dimen"),type="image",loader = "ectracker_popup.gif"),
                                            
                                            
                                  )
                                  
                                ),
                       tabPanel(tags$b("2. scData Processing"),
                                
                                
                                actionButton("action1","One-Click Analysis*",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                                br(),br(),
                                h6("Object"),
                             withLoader(verbatimTextOutput("seu_object"),type="image",loader = "ectracker_popup.gif"),br(),
                                actionButton("EX1","Explanation",icon = icon("book-open")),
                                withLoader(verbatimTextOutput("ex_object"),type="image",loader = "ectracker_popup.gif"),
                       
                       br(),br(),
                    
                                h6("Variable Feature Plot"),
                                
                               withLoader( plotOutput("p1"),type="image",loader = "ectracker_popup.gif"),
                                
                                br(), downloadButton(outputId="dndPlot","Download Plot"),    
                     actionButton("EX2","Explanation", icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_variable"),type="image",loader = "ectracker_popup.gif"),
                     
                    
                    
                    br(),br(),
                      h6("Variable Features"),
                      withLoader(plotOutput("plot_vf"), type="image", loader = "ectracker_popup.gif"),
                     br(), downloadButton(outputId="vfPlot","Download Plot"),
                     actionButton("EX4","Explanation",icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_vfeature"),type="image",loader = "ectracker_popup.gif"),
                    
                     
                     
                      
                      br(),br(),
                      h6("PCA"),
                      withLoader(plotOutput("pca"),type="image",loader = "ectracker_popup.gif"),
                      br(),downloadButton(outputId="pcaPlot","Download Plot"),
                      actionButton("EX6","Explanation", icon = icon("book-open")),
                       withLoader(verbatimTextOutput("ex_pca"),type="image",loader = "ectracker_popup.gif"),
                    
                      
                     br(),br(),
                    
                      h6("Heat Map"),
                      withLoader(plotOutput("heatmap"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                      
                      downloadButton(outputId="heatmapPlot","Download Plot"), actionButton("EX7","Explanation", icon = icon("book-open")),
                      
                      withLoader(verbatimTextOutput("ex_heatmap"),type="image",loader = "ectracker_popup.gif"),
                      
                      
                     br(), 
                     br(),
                     h6("JackStraw Plot"),
                      
                     withLoader(plotOutput("jackstraw"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                    
                     downloadButton(outputId="jackstrawPlot","Download Plot"),actionButton("EX8","Explanation", icon = icon("book-open")),
                     
                     withLoader(verbatimTextOutput("ex_jackstraw"),type="image",loader = "ectracker_popup.gif"),
                      
                     br(),  
                     br(),
                     h6("Elbow Plot"),
                      withLoader(plotOutput("elbow"),type="image",loader = "ectracker_popup.gif"),
                     br(),
                    downloadButton(outputId="elbowPlot","Download Plot"),actionButton("EX9","Explanation", icon = icon("book-open")),
                    
                    withLoader(verbatimTextOutput("ex_elbow"),type="image",loader = "ectracker_popup.gif"),
                    br(),br(),

                    
                      
                      h6("UMAP"),
                      withLoader(plotOutput("umap"),type="image",loader = "ectracker_popup.gif"),
                      br(),
                      downloadButton(outputId="umapPlot","Download Plot"),
                      actionButton("EX10","Explanation", icon = icon("book-open")),
                      
                      withLoader(verbatimTextOutput("ex_umap"),type="image",loader = "ectracker_popup.gif"),
                    actionButton("toTop", "Top", icon = icon("arrow-alt-circle-up")),
                    
                    ),
                      
                    
                    tabPanel(tags$b("3. DEG Analysis"),
                             sliderInput("log_fc", "Select Fold Change Threshold (log)",
                                         min = 0, max = 20,step=0.25,
                                         value = 1),     
                     downloadLink('downloadData', 'Download'),
                      br(),br(),
                     h6("Markers Table"),
                    withLoader(tableOutput("markers_display"),type="image",loader = "ectracker_popup.gif"),
                     useShinyjs(),
                     br(), 
                    fileInput("meta","Upload Meta File"),
                    checkboxInput("testme3", " Transdifferentiation (Meta file)", 
                                  value = FALSE),br(),
                    actionButton("m1","Table with Meta File "), 
                    withLoader(verbatimTextOutput("merge_meta"),type="image",loader = "ectracker_popup.gif"),
                    downloadLink('downloadData_meta', 'Download'),
                     br(),
                   br(),
                    h6("Generate Cluster Table"),
                      withLoader(verbatimTextOutput("cluster_ident"),type="image",loader = "ectracker_popup.gif"),
                     downloadLink('downloadData_cluster_ident', 'Download'),br(),br(),
                     extendShinyjs(text = jscode,functions = c("toTop")),
                    actionButton("toTop1", "Top", icon = icon("arrow-alt-circle-up")),
                    
                     
                     
                    ),
                   tabPanel(tags$b("4. CellEnrich "),
                            
                            tabsetPanel(type = "tabs",
                                        tabPanel(" AUCell (Method I)",
                                                 br(),br(),br(),
                                                  
                                                 br(),
                                                 selectInput("aucell","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Embryonic_stem_cell","Adult b_cell","Adult b_cell_plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal_cell","Adult CB CD34+23","Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 
                                                 br(), h6("Cell Ranking"), 
                                                 withLoader(plotOutput("cell_rank"),type="image",loader = "ectracker_popup.gif"),
                                                 br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_cellrank"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 br(), br(),   h6(" Cell Assignment"),
                                                 withLoader(plotOutput("au_cell"),type="image",loader = "ectracker_popup.gif"),br(),
                                                 actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_assign"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 br(), br(),  h6("UMAP"),
                                                 withLoader(plotOutput("tsne"),type="image",loader = "ectracker_popup.gif"),
                                                 br(),  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_uma"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                              br(),   br(),  h6("Plot for Signatures"),
                                                 withLoader(plotOutput("c_tsne"),type="image",loader = "ectracker_popup.gif"),
                                                 br(),
                                                  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_umappp"),type="image",loader = "ectracker_popup.gif"),
                                                downloadButton(outputId="cell_assign","Download All Plots"),                                                 
                                                downloadLink(outputId="cell_assign_score","Download Cell Assignment Table"),
                                                
                                                
                                                 br(), br(),  actionButton("toTop3", "Top", icon = icon("arrow-alt-circle-up")),
                                                 
                                                 
                                        ),
                                        
                                        tabPanel(" Stouffer Score (Method II)",
                                                 br(),br(),br(),br(),
                                                 selectInput("st_sc","Select Signature",choices=c("Devloping heart atrial cardiomyocyte","Devloping heart ventricular cardiomyocyte","Embryonic_stem_cell","Ectodermal cell differentiation","Ectodermal development","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Myofibroblast differentiation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult b_cell","Adult b_cell_plasmocyte","Adult basal_cell","Adult CB CD34+23","Adult chondrocyte","Adult chemosensory ","Adult dendritic_cell ","Adult endothelial_cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte_progenitor","Adult enterocyte","Adult epithelial cell (intermediated)","Adult epithelial_cell","Adult erythroid_cell ","Adult erythroid progenitor cell (RP high)","Adult fasciculata_cell ","Fetal B cells","Fetal Basophil_Mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal Collecting duct lineage","Fetal Connecting tubule lineage","Fetal Distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal Enteroendocrine cells_intestine_pancreas","Fetal Enteroendocrine cells_stomach","Fetal Erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ Microglia","Fetal ILC 3","Fetal Islet beta cells","Fetal Islet delta cells","Fetal LEC","Fetal Loop of Henle lineage","Fetal Macrophages","Fetal Meg progenitors","Fetal Meg","Fetal Megakaryoblasts","Fetal Microglia","Fetal Nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal Perivascular macrophages","Fetal Phagocytic macrophages","Fetal Plasma cells","Fetal Podocyte lineage","Fetal Proximal tubule lineage","Fetal PTPRC+ Microglia","Fetal Pulmonary neuroendocrine cells","Fetal Renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ Microglia","Fetal TRAF1+ APCs","Fetal Ureteric tip","Fetal Ureteric trunk","Fetal Antigen-presenting macrophages","Fetal Endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal fibroblast","Fetal skeletal muscle cell","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet_cell ","Adult hepatocyte_Endodermal cell","Adult hESC87","Adult immature sertoli cell (Pre-Sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 Macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal_cell","Adult smooth_muscle_cell","Adult stratified_epithelial_cell","Adult T_cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 br(),
                                                 
                                                 h6("Stouffer Plot"),
                                                 withLoader(plotOutput("stouffer_sc_plot"),type="image",loader = "ectracker_popup.gif"),
                                                 br(), downloadButton(outputId="cell_stoudownload","Download Plot"),
                                                 
                                                 actionButton("ex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 br(),   br(),h6("One Sided Wilcoxon P-Value"),
                                                 withLoader(uiOutput("Stouffer_sc"),type="image",loader = "ectracker_popup.gif"), br(),
                                                downloadLink('stouffer_table_cell', 'Download Table'),br(),br(),
                                                 
                                                 
                                                 
                                                 
                                        ) 
                                        
                                        
                                        
                                        
                            )),
                   
                   tabPanel( tags$b("5. TissueEnrich"),
                             tabsetPanel(type = "tabs",
                                         
                                         
                                         
                                         tabPanel("Hypergeometric method (Method I)",
                                                  br(),
                                                  sliderInput("obs", "Select Cluster Number ",
                                                              min = 0, max = 20,
                                                              value = 0),
                                                   downloadLink('downloadData_cluster_name', 'Download Cluster Table'),br(),br(),
                                                   
                                                  br(),br(),   numericInput("dataset", "Enter Dataset Number: 1: HPA ; 2: GTEx", 1, min = 1, max = 2), br(),
                                                  
                                                  h6("Generate Tissue Enrich Table"),
                                                  
                                                  withLoader(verbatimTextOutput("tissue_detail"),type="image",loader = "ectracker_popup.gif"),
                                                  
                                                  downloadLink('downloadtissuedata ', 'Download'),
                                                  
                                                  br(),  
                                                  br(),
                                                  
                                                 h6("Fold Change Plot"),
                                                  
                                                  withLoader(plotOutput("plot_fold"),type="image",loader = "ectracker_popup.gif"), 
                                                  br(),
                                                  downloadButton(outputId="tissue_fold1","Download Plot"),
                                                  
                                                  
                                                  br(),br(),
                                                  
                                                  
                                                  h6(" Log10(PValue)Plot"),
                                                  
                                                  withLoader(plotOutput("plot_log"),type="image",loader = "ectracker_popup.gif"),
                                                  br(),   downloadButton(outputId="tissue_log1","Download Plot"),
                                                  
                                                  
                                                  br(),br(),
                                                  
                                                  
                                                  h6("Tissue Specific Visualization"),
                                                  br(), uiOutput("t_or"),
                                                  uiOutput("tt_or"),br(),br(),
                                                  br(), plotOutput("tor_umap"),br(),
                                                  downloadButton(outputId="torPlot","Download Plot"),
                                                  actionButton("EX101","Explanation", icon = icon("book-open")),br(),
                                                  br(),br(),    
                                                  h6("Genomic Location"),br(),br(),
                                                  withLoader(plotOutput("de_karyotype"),type="image",loader = "ectracker_popup.gif"),
                                                  br(), 
                                                 downloadButton(outputId="karyoplott","Download Plot"),
                                                 actionButton("EX2","Explanation", icon = icon("book-open")),
                                                 
                                                 withLoader(verbatimTextOutput("ex_karyo"),type="image",loader = "ectracker_popup.gif"),
                                                 
                                                 
                                                 br(),
                                                 
                                                 
                                                 
                                                 
                                                 actionButton("toTop4", "Top", icon = icon("arrow-alt-circle-up")),
                                       
                                                   ),
                                         
                                         
                                         tabPanel("Stouffer Score (Method II)",
                                                  
                                                  tabsetPanel( type= "tabs", 
                                                               
                                                               tabPanel("HPA Data",
                                                                        
                                                                        br(),br(),

                                                                        selectInput("hpa11","Select Tissue",choices=c("Lymph Node",	"Tonsil",	"Appendix",	"Spleen",	"Bone Marrow"	,"Esophagus"	,"Skin",	"Colon",	"Rectum",	"Duodenum",	"Small Intestine",	"Stomach",	"Adipose Tissue",	"Lung",	"Placenta",	"Gallbladder",	"Urinary Bladder","Endometrium",	"Smooth Muscle","Fallopian Tube",	"Thyroid Gland",	"Ovary","Prostate",	"Kidney",	"Adrenal Gland",	"Brain",	"Salivary Gland","Pancreas",	"Skeletal Muscle",	"Liver",	"Heart Muscle")),
                                                                        
                                                                        
                                                                        h6("Stouffer Score "),
                                                                        withLoader(plotOutput("hpa_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                        br(),  
                                                                      downloadButton(outputId="hpaplot","Download Plot"),
                                                                      actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                      withLoader(verbatimTextOutput("ex_stouffer_hpa"),type="image",loader = "ectracker_popup.gif"),
                                                                      
                                                                      
                                                                      br(),br(),
                                                                        h6("One Sided Wilcoxon P-Value"),
                                                                        withLoader(uiOutput("hpa_set2"),type="image",loader = "ectracker_popup.gif"),
                                                                      downloadLink('stouffer_table_hpa', 'Download Table'),br(),br(),
                                                                        actionButton("toTop5", "Top", icon = icon("arrow-alt-circle-up")),
                                                                        
                                                                        
                                                               ),
                                                               
                                                               
                                                               br(), tabPanel("GTEx Data",br(),br(),
                                                                              selectInput("gtex11","Select Tissue",choices=c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Blood Vessel","Uterus","Cervix/Uterine","Pituitary","Vagina","Pancreas")),
                                                                              
                                                                              br(),h6("Stouffer Score"),
                                                                              withLoader(plotOutput("gtex_stouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                              br(), 
                                                                              downloadButton(outputId="gtexplot","Download Plot"),
                                                                              actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                              withLoader(verbatimTextOutput("ex_stouffer_gtex"),type="image",loader = "ectracker_popup.gif"),
                                                                              
                                                                              br(),br(),
                                                                              
                                                                              h6("One Sided Wilcoxon P-Value"),
                                                                              withLoader(uiOutput("gtex_pstouffer"),type="image",loader = "ectracker_popup.gif"),
                                                                              downloadLink('stouffer_table_gtex', 'Download Table'),br(),br(),
                                                                              
                                                                              
                                                                              actionButton("toTop100", "Top", icon = icon("arrow-alt-circle-up")),
                                                               ) )             
                                         )
                                         
                                         
                             )),
                   
                   
                   tabPanel(tags$b("6. Gene Regulatory Network"),
                            tabsetPanel(
                              tabPanel("Cell Cluster",
                                       br(),br(),
                                       sliderInput("grn_obs", "Select Cluster Number ",
                                                   min = 0, max = 20,
                                                   value = 0),
                                       h6("Generate Specific Cluster Table"),br(),
                                       
                                       
                                       withLoader(verbatimTextOutput("grn_clust"),type="image",loader = "ectracker_popup.gif"),  
                                       
                              ),
                              
                              
                              tabPanel("CellEnrich",
                                       
                                       br(),br(),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Transcription Factors and Target Genes)"),
                                      withLoader(simpleNetworkOutput("cell_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Target Genes and Tissue Signatures)"),
                                      withLoader(simpleNetworkOutput("cell_grn"),type="image",loader = "ectracker_popup.gif"),
                                      

                                       actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                       withLoader(verbatimTextOutput("ex_cellgrn"),type="image",loader = "ectracker_popup.gif"),
                                       
                                       
                                       br(),br(),
                                       h6("Gene Regulatory Table"),
                                       withLoader(uiOutput("cell_grn_table"),type="image",loader = "ectracker_popup.gif"),br(),
                                       actionButton("ex_cellgrn","Explanation", icon = icon("book-open")), downloadLink('download_cell_gt', 'Download Table'),br(),br(),

                                       withLoader(verbatimTextOutput("ex_cellgrn_table"),type="image",loader = "ectracker_popup.gif"),
                                       br(), actionButton("toTop6", "Top", icon = icon("arrow-alt-circle-up")),
                                       
                              ),tabPanel("TissueEnrich",
                                         
                                         br(),br(),

                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("hpa_grn"),type="image",loader = "ectracker_popup.gif"),
                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("hpa_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         br(),
                                         
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         h6("Gene Regulatory Table (HPA)"),
                                         withLoader(uiOutput("hpa_grn_table"),type="image",loader = "ectracker_popup.gif"),
                                         br(),
                                         downloadLink('download_hpa_gt', 'Download Table'),br(),br(),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn_table"),type="image",loader = "ectracker_popup.gif"),br(),br(),
                                         
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("gtex_grn_1"),type="image",loader = "ectracker_popup.gif"),
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("gtex_grn"),type="image",loader = "ectracker_popup.gif"),
                                         br(),
                                         
                                         
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_gtexgrn"),type="image",loader = "ectracker_popup.gif"),
                                         
                                         br(), h6("Gene Regulatory Table (GTEx)"),
                                         withLoader(uiOutput("gtex_grn_table"),type="image",loader = "ectracker_popup.gif"),
                                         br(),br(),br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),downloadLink('download_gtex_gt', 'Download Table'),
                                         withLoader(verbatimTextOutput("ex_gtexgrn_table"),type="image",loader = "ectracker_popup.gif"),
                                         actionButton("toTop7", "Top", icon = icon("arrow-alt-circle-up")),
                              )
                              
                              
                            )
                            
                   )               
                   
                   
                   
                   
                   
                    )
                      
                      
             ),
             
            
             
             
             
         
             
            
             
            
             tabPanel("Tutorial", icon = icon("angle-double-right"),
                     h6( "Check Out the Video Tutorial here",align="center"),
                    tags$iframe(width="1500", height="650", src="atlast_final.mp4", frameborder="0", allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"),
                   h6("Check out the tutorial page link here:"), tags$a(href="final_documentation.html", "Tutorial.html", style = "color:yellow"),
                   br(),br(),br(),br(),br(),

                      ),
             
            
             
             
             
             tabPanel("Contact Us", icon = icon("angle-double-right"),
   tabsetPanel( type = "tabs",
               
                         
                tabPanel("Developer",
                         HTML('<left><img src="vishakha.jpg", width="20%", height= "20%" ></left>'),
h6("Vishakha"), 
h6("Ph.D Scholar"),
h6("Email: vishakhag@iiitd.ac.in
"),
               ), 

tabPanel("Developer",
         HTML('<left><img src="siddhant.jpg", width="20%", height= "20%" ></left>'),
         h6("Siddhant"), 
         h6("M.tech student"),
         h6("Email: siddhant18241@iiitd.ac.in
")
          ),
tabPanel("Web Designer",
         HTML('<left><img src="aayushi.jpg", width="20%", height= "20%" ></left>'),
         h6("Aayushi Mittal"), 
         h6("Ph.D Scholar"),
         h6("Email: aayushim@iiitd.ac.in
")   
),
tabPanel("Testing",
         HTML('<left><img src="sanjay.jpg", width="20%", height= "20%" ></left>'),
         
         h6("Sanjay kumar Mohanty"), 
         h6("Ph.D Scholar"),
         h6("Email: sanjaym@iiitd.ac.in
")   ),
tabPanel("Collaborator",
         HTML('<left><img src="sengupta.PNG", width="20%", height= "20%" ></left>'),
         
       br(),br(),  tags$a(href="https://www.debarka.com/", "Sengupta Lab", style = "color:yellow")
         
         
         
         
         )





))

)
)
)
