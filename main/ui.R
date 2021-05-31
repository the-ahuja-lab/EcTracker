

jscode <- "shinyjs.toTop = function() {window.scrollTo(0, 0);}" ### for to top button


shinyUI(fluidPage(
  theme=shinytheme("cyborg"),
  tags$style(HTML("
        .tabs-above > .nav > li[class=active] > a {
           background-color: #000;
           color: #FFFF;
        }")),

 useShinyalert(),

  navbarPage("EcTracker",
                   tabPanel("Home" ,icon = icon("home"),

                           HTML('<center><img src="main.gif", width="58%", height= "58%" ></center>'),



                          column(8, align="center", offset = 2 ,print("Developed by: "), tags$u("Vishakha Gautam"), print("&" ), tags$u("Siddhant Kalra"),
                                 br(),print("For more details:"),
                                 tags$head(tags$style(HTML("a {color: blue}"))),
                                 tags$a(href="https://ahuja-lab.in/", "The Ahuja Lab", style = "color:yellow"))
                          
                          
                   ),

         
             tabPanel("Analysis", icon = icon("angle-double-right"),

                     navlistPanel(
                     tabPanel(tags$b("1. Upload") ,

                                  mainPanel(
                                    h5("A."),    h6("Sample Data"),
                                    checkboxInput("testme", "Sample Data ",
                                                  value = FALSE),
                                    checkboxInput("testme3", "Sample Metadata File",
                                                  value = FALSE),
                                    h6("..............................................................................................................................................................................................................................................................................................................................."),br(),
                       
                                    h5("B."),   h6(" Upload Expression Matrix"),
                                    column(width=4,  fileInput("file","Upload Expression Matrix (Maximum Size Limit is 10 GB)")),
                                   
                                     tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73122;", "GSE73122", style = "color:yellow"),br(),
                                    
                                    h6("..............................................................................................................................................................................................................................................................................................................................."),br(),
                                    
                                    h5("C.") ,    h6("Upload 10x Genomics Data Files"),
                                    
                                    fluidRow(column(width=4, fileInput("File100"," Upload barcodes.tsv File")),column(width=4,fileInput("File101","Upload genes.tsv File")),  column(width=4,fileInput("File102","Upload matrix.mtx File"))),
                                    
                                    downloadLink('barcode', 'Download barcode.tsv,', class = "butt",style = "color:yellow"),
                                    
                                    
                                    downloadLink('genes', 'genes.tsv,', class = "butt",style = "color:yellow"),
                                   
                                    downloadLink('matrix', 'matrix.mtx', class = "butt",style = "color:yellow"),
                                   
                                    h6("..............................................................................................................................................................................................................................................................................................................................."),br(),
                                    
                                    

                                  h5("D."),  h6("Upload Metadata File"),    
                                    
                           
                                    fluidRow( column(width=4, fileInput("meta","Upload Metadata File "))),
                                    
                                   
                                    br(),
                                    h6("..............................................................................................................................................................................................................................................................................................................................."),br(),
                                    
                                             h6("Information About Data and Dimensions"),
                                           withLoader(tableOutput("dat"),type="image",loader = "new_loader.gif"),
                                           
                                           withLoader(verbatimTextOutput("dimen"),type="image",loader = "new_loader.gif"),
                                    br(),br(),
                                    
                                    actionButton("toTop2025", "Top", icon = icon("arrow-alt-circle-up")),            

                                  )

                                ),

                       tabPanel(tags$b("2. scData Processing"),

tabsetPanel(type = "tabs",

            
            
            tabPanel("Method 1 (faster & abstract)",
                     br(),
                     checkboxInput("ERCC2", "Exclude Spike-Ins/Heterogenous Genes",
                                   value = FALSE),
                  
                     checkboxInput("low_seq_depth1", "Filter Cells with Low Sequencing Depth (Require Mapped_reads Information in metadata file)",
                                   value = FALSE),
                     sliderInput("min_seqdepth1", "Minimum Sequencing Depth of Cells",
                                 min = 20000, max = 100000,step=1000,
                                 value = 50000),
                     sliderInput("hvg", "Number of Features to Select as Top Variable Features",
                                 min = 500, max = 5000,step=1,
                                 value = 5000),
                     br(),
                     actionButton("SCE1","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                     br(),br(),
                      h6("QC-Plot "),
                     
                     withLoader(plotOutput("SCE_plot"),type="image",loader = "Asset 1.jpg"),
                     br(), downloadButton(outputId="qc_scater_download","Download Plot"),
                     actionButton("q1","Explanation",icon = icon("book-open")),
                     withLoader(verbatimTextOutput("q1_explanation"),type="image",loader = "Asset 1.jpg"),

                     br(),br(),
                     h6("Cells after Quality Control"),
                    withLoader(plotOutput("qc_scater"),type="image",loader = "Asset 1.jpg"),
                  br(),  downloadButton(outputId="qc_scater_download1","Download Plot"),
                    actionButton("q2","Explanation",icon = icon("book-open")),
                    withLoader(verbatimTextOutput("q2_explanation"),type="image",loader = "Asset 1.jpg"),
                    br(),br(),
                  h6("Cells after Quality Control"),
                  withLoader(verbatimTextOutput("cells_qc"),type="image",loader = "Asset 1.jpg"),
                  
                    br(),br(),
                    h6("t-SNE Plot"),
                    withLoader(plotOutput("SCE_tsneplot"),type="image",loader = "Asset 1.jpg"),
                    br(), downloadButton(outputId="drop_clust_scatter","Download Plot"),
                    actionButton("DEX1","Explanation",icon = icon("book-open")),
                    withLoader(verbatimTextOutput("dex_scatter"),type="image",loader = "Asset 1.jpg"),br(),br(),
                    
                     h6("Clustering"),
                     withLoader( plotOutput("SCE_plot_clustering"),type="image",loader = "Asset 1.jpg"),
                     br(), downloadButton(outputId="sce_clus_d","Download Plot"),
                     actionButton("Dsce","Explanation",icon = icon("book-open")),
                     withLoader(verbatimTextOutput("dsce_scatter"),type="image",loader = "Asset 1.jpg"),
                     br(),br(),
                     h6("Silhouette Plot"),
                     withLoader( plotOutput("SCE_plot_silhouette"),type="image",loader = "Asset 1.jpg"),
                     br(), downloadButton(outputId="sce_silhouette_d1","Download Plot"),
                     actionButton("Dsce1","Explanation",icon = icon("book-open")),
                     withLoader(verbatimTextOutput("dsce_silhouette"),type="image",loader = "Asset 1.jpg"),
                    br(),br(),
                    
                    actionButton("toTop2021", "Top", icon = icon("arrow-alt-circle-up")),

                     
            ),
            

            tabPanel("Method 2 (slower & comprehensive)",
                     

                     br(),
                                
                     sliderInput("min_cell", "Minimum Number of Cells",
                                 min = 3, max = 100,step=1,
                                 value = 3),
                     sliderInput("min_features", "Minimum Number of Features",
                                 min = 200, max = 1000,step=1,
                                 value = 200),
                     checkboxInput("ERCC", "Exclude Spike-Ins/Heterogenous Genes",
                                   value = FALSE),
                     checkboxInput("low_seq_depth", "Filter Cells with Low Sequencing Depth (Require Mapped_reads Information in metadata file)",
                                   value = FALSE),
                     sliderInput("min_seqdepth", "Minimum Sequencing Depth of Cells",
                                 min = 20000, max = 100000,step=1000,
                                 value = 50000),
                     br(),br(),
                     actionButton("seurat1","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                     br(),br(),
                                h6("Object"),
                             withLoader(verbatimTextOutput("seu_object"),type="image",loader = "Asset 1.jpg"),br(),
                                actionButton("EX1","Explanation",icon = icon("book-open")),
                                withLoader(verbatimTextOutput("ex_object"),type="image",loader = "Asset 1.jpg"),

                       br(),br(),

                                h6("QC Plot"),

                               withLoader( plotOutput("p1"),type="image",loader = "Asset 1.jpg"),

                                br(), downloadButton(outputId="dndPlot","Download Plot"),
                     actionButton("EX2","Explanation", icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_variable"),type="image",loader = "Asset 1.jpg"),
                     br(),br(),
                     
                     h6("Scatter Plot of Feature-Feature Relationship"),
                     
                     withLoader( plotOutput("featurescatter_plot"),type="image",loader = "Asset 1.jpg"),
                     br(),downloadButton(outputId="featurescatterPlot","Download Plot"),
                     actionButton("featurescatter_explanation1","Explanation", icon = icon("book-open")),
                     withLoader(verbatimTextOutput("featurescatter_explanation"),type="image",loader = "Asset 1.jpg"),
                     br(),br(),
                     
                    
                     sliderInput("nFeature_rna", "nFeature_RNA Minimum Value",
                                 min = 200, max = 30000,step=1,
                                 value = 200),
                     sliderInput("count_rna", "nFeature_RNA Maximum Value",
                                 min = 1000, max = 30000,step=1,
                                 value = 10000),
                     sliderInput("mito_rna", "percent.mt Maximum Value",
                                 min = 0, max = 100,step=1,
                                 value = 5),
                     sliderInput("ercc_rna", "percent.ERCC Maximum Value",
                                 min = 0, max = 100,step=1,
                                 value = 5),
                     br(),br(),
                     
                     h6("Cells after Quality Control"),
                     withLoader(verbatimTextOutput("qc_seurat"),type="image",loader = "Asset 1.jpg"),
                     
                     br(),br(),
                     
                     
                     sliderInput("find_var_feature", "Number of Features to Select as Top Variable Features",
                                 min = 500, max = 5000,step=1,
                                 value = 2000),
                      h6("Variable Features"),
                      withLoader(plotOutput("plot_vf"), type="image", loader = "Asset 1.jpg"),
                     br(), downloadButton(outputId="vfPlot","Download Plot"),
                     actionButton("EX4","Explanation",icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_vfeature"),type="image",loader = "Asset 1.jpg"),
                    br(),br(),
                    checkboxInput("vars_regress", "Regress out Variables",
              value = FALSE),



                      br(),br(),
                      h6("PCA"),
                      withLoader(plotOutput("pca"),type="image",loader = "Asset 1.jpg"),
                      br(),downloadButton(outputId="pcaPlot","Download Plot"),
                      actionButton("EX6","Explanation", icon = icon("book-open")),
                       withLoader(verbatimTextOutput("ex_pca"),type="image",loader = "Asset 1.jpg"),


                     br(),br(),

                      h6("Heat Map"),
                      withLoader(plotOutput("heatmap"),type="image",loader = "Asset 1.jpg"),
                      br(),

                      downloadButton(outputId="heatmapPlot","Download Plot"), actionButton("EX7","Explanation", icon = icon("book-open")),

                      withLoader(verbatimTextOutput("ex_heatmap"),type="image",loader = "Asset 1.jpg"),


                     br(),
                     br(),
                     h6("JackStraw Plot"),

                     withLoader(plotOutput("jackstraw"),type="image",loader = "Asset 1.jpg"),
                      br(),

                     downloadButton(outputId="jackstrawPlot","Download Plot"),actionButton("EX8","Explanation", icon = icon("book-open")),

                     withLoader(verbatimTextOutput("ex_jackstraw"),type="image",loader = "Asset 1.jpg"),

                     br(),
                     br(),
                     h6("Elbow Plot"),
                      withLoader(plotOutput("elbow"),type="image",loader = "Asset 1.jpg"),
                     br(),
                    downloadButton(outputId="elbowPlot","Download Plot"),
                    actionButton("EX9","Explanation", icon = icon("book-open")),

                    withLoader(verbatimTextOutput("ex_elbow"),type="image",loader = "Asset 1.jpg"),
                    br(),br(),



                      h6("UMAP"),
                      withLoader(plotOutput("umap"),type="image",loader = "Asset 1.jpg"),
                      br(),
                      downloadButton(outputId="umapPlot","Download Plot"),
                      actionButton("EX10","Explanation", icon = icon("book-open")),
                    withLoader(verbatimTextOutput("ex_umap"),type="image",loader = "Asset 1.jpg"),
                    
                    br(),br(),
                    
                    actionButton("toTop", "Top", icon = icon("arrow-alt-circle-up")),

                    
                    
                    
                    
                    
                    )
            
            )


),
tabPanel(tags$b("3. Batch Correction (optional)"),
         
         tabsetPanel(type = "tabs",
                     
                     
                     tabPanel("Method 1 (faster & abstract)",
                              
                             br(), h6("For this step, please upload the metadata file in Upload Section"),br(),br(),
                              
                              actionButton("Batch_Correction_method_dropclust1","Merge Metadata File"), br(),br(),
                              withLoader( uiOutput("meta_info_dropclust"),type="image",loader = "Asset 1.jpg"),
                              br(),br(),
                              actionButton("dropclust_batch1","Batch Correction",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
                              
                              withLoader( plotOutput("sce_batch"),type="image",loader = "Asset 1.jpg"),
                             br(), downloadButton(outputId="drop_clust_batch_Plot","Download Plot"),
                              actionButton("DEX_batch2","Explanation", icon = icon("book-open")),
                              withLoader(verbatimTextOutput("dex_scatter1"),type="image",loader = "Asset 1.jpg"),
                              br(),br(),
                              h6("Corrected Expression for Gene "),br(), 
                              sliderInput("sce_gene", "Select Gene Number",
                                          min = 0, max = 30000,step=1,
                                          value = 1),
                              withLoader( plotOutput("drop_clust_markers_batch"),type="image",loader = "Asset 1.jpg"),                               
                             br(), downloadButton('drop_clust_batch_Plot_hist', 'Download'),
                              actionButton("DEX_batch3","Explanation", icon = icon("book-open")),
                              withLoader(verbatimTextOutput("dex_scatter2"),type="image",loader = "Asset 1.jpg"),
                             br(),br(),
                             
                             actionButton("toTop1000101", "Top", icon = icon("arrow-alt-circle-up")),
                             ),
                     
                     tabPanel("Method 2 (slower & comprehensive)",
                              
                              
                                           
                                   tabPanel("Harmony",
                                           
                                            br(),
                                            h6("For this step, please upload the metadata file in Upload Section"),br(),br(),
                                            actionButton("Batch_Correction_method1","Merge Metadata File"),br(),br(),
                                            withLoader( uiOutput("meta_info"),type="image",loader = "Asset 1.jpg"),
                                            
                                            br(),br(),
                                            actionButton("Batch_Correction_method","Batch Correction",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")), br(),br(),
                                            h6("Harmony Run"),
                                            withLoader(plotOutput("harm1"),type="image",loader = "Asset 1.jpg"),
                                            br(),
                                            downloadButton(outputId="harmony_plot1","Download Plot"),
                                            actionButton("hEX2","Explanation", icon = icon("book-open")),
                                            withLoader(verbatimTextOutput("hex_harmplot"),type="image",loader = "Asset 1.jpg"),
                                            
                                            br(),br(),
                                            h6("Embeddings"),
                                            withLoader(verbatimTextOutput("harm2"),type="image",loader = "Asset 1.jpg"),
                                            
                                            downloadLink('harm_tab', 'Download', class = "butt",style = "color:yellow"),
                                            
                                            br(),br(),br(),
                                            
                                            h6("Dimensionality Reduction and Variable Feature plot"),
                                            withLoader(plotOutput("harm3"),type="image",loader = "Asset 1.jpg"),
                                            br(),
                                            downloadButton(outputId="harmony_plot2","Download Plot"),
                                            actionButton("hEX1","Explanation", icon = icon("book-open")),
                                            withLoader(verbatimTextOutput("hex_umap"),type="image",loader = "Asset 1.jpg"),
                                            
                                            br(),br(),
                                            
                                            h6("UMAP"),
                                            withLoader(plotOutput("harm4"),type="image",loader = "Asset 1.jpg"),
                                            br(),
                                            downloadButton(outputId="harmony_plot3","Download Plot"),
                                            actionButton("hEX9","Explanation", icon = icon("book-open")),
                                            withLoader(verbatimTextOutput("hex_pca"),type="image",loader = "Asset 1.jpg"),
                                            br(),br(),
                                            
                                            actionButton("toTop10001010", "Top", icon = icon("arrow-alt-circle-up")),
                                             )
                                   
                            
                              
                               )
                     
)),





                    tabPanel(tags$b("4. DGE Analysis"),
                             
                             tabsetPanel(type = "tabs",
                            
                                         tabPanel    ("Batch Uncorrected",
                             br(),
                             
                             tabsetPanel(type = "tabs",
                                         

                                         tabPanel("Method 1 (faster & abstract)",
                               br(),br(),
                               actionButton("scater_de_action1","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),br(),br(),br(),
                               h6("Differential Gene Expression"),br(),
                               selectInput("test_use_sce","Select Test for DE Genes",choice=c("Pairwise Wilcoxon Rank Sum Test","Pairwise t-Test","Pairwise Binomial Test")),
                              
                              withLoader( verbatimTextOutput("SCE_markers"),type="image",loader = "Asset 1.jpg"),
                               
                            br(),br(),

                               h6("Display Markers"),
                              
                            
                               withLoader(verbatimTextOutput("SCE_marker_display"),type="image",loader = "Asset 1.jpg"),
                               downloadLink('drop_clust_de', 'Download', class = "butt",style = "color:yellow"),
                               br(),br(),
                               h6("Feature Plot"),
                            withLoader( uiOutput("de_f"),type="image",loader = "Asset 1.jpg"),

                              withLoader(plotOutput("feature_plot_result_sce"),type="image",loader = "Asset 1.jpg"),
                            br(),br(),
                            
                            actionButton("toTop2022", "Top", icon = icon("arrow-alt-circle-up")),
                    ),
                  
                    tabPanel("Method 2 (slower & comprehensive)", 
                             br(), sliderInput("log_fc", "Select Fold Change Threshold (log)",
                                                    min = 0, max = 5,step=0.25,
                                                    value = 1),
                             selectInput("test_use","Select Test for DE Genes",choice=c("Wilcoxon","Bimod (Likelihood-ratio test for single cell gene expression)","ROC analysis","Student's t-test","Negative binomial generalized linear model","Poisson generalized linear model","Logistic regression ","MAST","DESeq2")),
                             br(),
                             actionButton("seurat_de_action1","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),br(),
                             br(),
                             
                             h6("Markers Table"),
                             withLoader(tableOutput("markers_display"),type="image",loader = "Asset 1.jpg"),
                             downloadLink('downloadData', 'Download', class = "butt",style = "color:yellow"),
                             
                             useShinyjs(),
                             br(),
                             
                            
                             br(),
                             h6("Feature Plot"),
                             withLoader( uiOutput("de_sf"),type="image",loader = "Asset 1.jpg"),
                             
                             #textInput("feature_plot","Enter Gene Name"),
                             withLoader(plotOutput("feature_plot_result"),type="image",loader = "Asset 1.jpg"),
                             br(),
                             downloadButton(outputId="torPlot","Download Plot"),
                             actionButton("fEX2","Explanation", icon = icon("book-open")),
                             
                             withLoader(verbatimTextOutput("ex_feature"),type="image",loader = "Asset 1.jpg"),
                             
                             br(), br(),
                             h6("Genomic Location"),
                             withLoader(plotOutput("de_karyotype"),type="image",loader = "Asset 1.jpg"),
                             br(),
                             downloadButton(outputId="karyoplott","Download Plot"),
                             actionButton("EX2","Explanation", icon = icon("book-open")),
                             
                             withLoader(verbatimTextOutput("ex_karyo"),type="image",loader = "Asset 1.jpg"),
                             
                             br(),br(),
                             h6("Generate Cluster Table"),
                             withLoader(verbatimTextOutput("cluster_ident"),type="image",loader = "Asset 1.jpg"),
                             downloadLink('downloadData_cluster_ident', 'Download', class = "butt",style = "color:yellow"),br(),br(),
                             extendShinyjs(text = jscode,functions = c("toTop")),
                             actionButton("toTop1", "Top", icon = icon("arrow-alt-circle-up")),
                             
                             
                             
                    )
                  
                  )),
                  tabPanel("Batch Corrected (optional)",
                           br(),
                           actionButton("mast1","Merge Metadata File"), br(),br(),
                           withLoader( uiOutput("meta_info_mast"),type="image",loader = "Asset 1.jpg"),
                           withLoader( uiOutput("meta_info_condition"),type="image",loader = "Asset 1.jpg"),
                           withLoader( uiOutput("meta_info_condition1"),type="image",loader = "Asset 1.jpg"),
                           br(),br(),
                           sliderInput("m1", "Frequency of Expressed Genes ",
                                       min = 0, max = 50,step= 0.1,
                                       value = 0.2),
                           sliderInput("m2", "Fold Change Threshold (log2) ",
                                       min = 0, max = 20,step= 0.25,
                                       value = 1),
                           sliderInput("m3", "False Discovery Rate ",
                                       min = 0, max = 1,step= 0.01,
                                       value = 0.05),br(),br(),
                          
                           actionButton("mast2","One-Click Analysis",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")), br(),br(),
                           h6("Markers Identification"),
                           
                           
                           withLoader(  tableOutput("condition_info"),type="image",loader = "Asset 1.jpg"),
                           downloadLink('downloadData_mast', 'Download', class = "butt",style = "color:yellow"),br(),br(),
                          
                          
                           )
                  
                  
                  )),

                   tabPanel(tags$b("5. CellEnrich "),

                            tabsetPanel(type = "tabs",
                                        tabPanel(" AUCell (Method I)",
                                                 br(),

                                                 
                                                 actionButton("action1","One-Click Enrichment",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("paper-plane")),
br(),br(),
                                                 selectInput("aucell","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),

                                                 br(), h6("Cell Ranking"),
                                                 withLoader(plotOutput("cell_rank"),type="image",loader = "Asset 1.jpg"),
                                                 br(), actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_cellrank"),type="image",loader = "Asset 1.jpg"),


                                                 br(), br(),   h6(" Cell Assignment"),
                                                 withLoader(plotOutput("au_cell_R"),type="image",loader = "Asset 1.jpg"),br(),
                                                 actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_assign"),type="image",loader = "Asset 1.jpg"),


                                                

                                              br(),   br(),  h6("UMAP for Signatures"),
                                                 withLoader(plotOutput("c_tsne"),type="image",loader = "Asset 1.jpg"),
                                                 br(),
                                                  actionButton("ex11","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_umappp"),type="image",loader = "Asset 1.jpg"),
                                              br(),br(),
                                                downloadButton(outputId="cell_assign","Download All Plots"),
                                                downloadLink(outputId="cell_assign_score","Download Cell Assignment Table", class = "butt",style = "color:yellow"),

                                                 br(), br(), br(),  actionButton("toTop3", "Top", icon = icon("arrow-alt-circle-up")),


                                        ),

                                        tabPanel(" Stouffer Based Enrichment Analysis (Method II)",
                                                 br(),
                                                 selectInput("st_sc","Select Signature",choices=c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")),
                                                 br(),

                                                 h6("Stouffer Plot"),
                                                 withLoader(plotOutput("stouffer_sc_plot"),type="image",loader = "Asset 1.jpg"),
                                                 br(), downloadButton(outputId="cell_stoudownload","Download Plot"),

                                                 actionButton("ex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_stouffer"),type="image",loader = "Asset 1.jpg"),

                                                 br(),   br(),h6("One Sided Wilcoxon P-Value"),
                                                 withLoader(uiOutput("Stouffer_sc"),type="image",loader = "Asset 1.jpg"), br(),
                                                downloadLink('stouffer_table_cell', 'Download Table', class = "butt",style = "color:yellow"),br(),br(),

                                                br(), br(),  actionButton("toTop2023", "Top", icon = icon("arrow-alt-circle-up")),


                                        )




                            )),

                   tabPanel( tags$b("6. TissueEnrich"),
                             tabsetPanel(type = "tabs",



                                         tabPanel("Hypergeometric method (Method I)",
                                                  br(),
                                                  sliderInput("obs", "Select Cluster Number ",
                                                              min = 0, max = 20,
                                                              value = 1),
                                                   downloadLink('downloadData_cluster_name', 'Download Cluster Table', class = "butt",style = "color:yellow"),br(),br(),

                                                  br(),br(),   numericInput("dataset", "Enter Dataset Number: 1: HPA ; 2: GTEx", 1, min = 1, max = 2), br(),

                                                  h6("Generate TissueEnrich Table"),

                                                  withLoader(verbatimTextOutput("tissue_detail"),type="image",loader = "Asset 1.jpg"),

                                                  # downloadLink('downloadtissuedata ', 'Download', class = "butt",style = "color:yellow"),

                                                 
                                                  br(),br(),

                                                 h6("Fold Change Plot"),

                                                  withLoader(plotOutput("plot_fold"),type="image",loader = "Asset 1.jpg"),
                                                  br(),
                                                  downloadButton(outputId="tissue_fold1","Download Plot"),
                                                 
                                                 actionButton("foldex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_fold_change"),type="image",loader = "Asset 1.jpg"),

                                                  br(),br(),


                                                  h6(" Log10 (P-Value) Plot"),

                                                  withLoader(plotOutput("plot_log"),type="image",loader = "Asset 1.jpg"),
                                                  br(),   downloadButton(outputId="tissue_log1","Download Plot"),actionButton("logdex12","Explanation", icon = icon("book-open")),
                                                 withLoader(verbatimTextOutput("ex_log_p_value"),type="image",loader = "Asset 1.jpg"),


                                              


                                                
                                                 br(),br(),
                                                 h6("Tissue-Specific Genes"),

                                                 uiOutput("t_or"),
                                                 uiOutput("tt_or"),br(),br(),
                                                 br(), 

                                                 actionButton("toTop4", "Top", icon = icon("arrow-alt-circle-up")),

                                                   ),


                                         tabPanel("Stouffer Based Enrichment Analysis (Method II)",

                                                  tabsetPanel( type= "tabs",

                                                               tabPanel("HPA Data",

                                                                        br(),br(),
                                                                        selectInput("hpa11","Select Tissue",choices=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate-HPA",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")),

                                                                        br(),

                                                                        h6("Stouffer Score "),
                                                                        withLoader(plotOutput("hpa_stouffer"),type="image",loader = "Asset 1.jpg"),
                                                                        br(),
                                                                      downloadButton(outputId="hpaplot","Download Plot"),
                                                                      actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                      withLoader(verbatimTextOutput("ex_stouffer_hpa"),type="image",loader = "Asset 1.jpg"),


                                                                      br(),br(),
                                                                        h6("One Sided Wilcoxon P-Value"),
                                                                        withLoader(uiOutput("hpa_set2"),type="image",loader = "Asset 1.jpg"),
                                                                      downloadLink('stouffer_table_hpa', 'Download Table', class = "butt",style = "color:yellow"),br(),br(),
                                                                        actionButton("toTop5", "Top", icon = icon("arrow-alt-circle-up")),


                                                               ),


                                                                tabPanel("GTEx Data",br(),br(),
                                                                              selectInput("gtex11","Select Tissue",choices=c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix-Uterine","Pituitary","Vagina","Pancreas")),

                                                                              br(),h6("Stouffer Score"),
                                                                              withLoader(plotOutput("gtex_stouffer"),type="image",loader = "Asset 1.jpg"),
                                                                              br(),
                                                                              downloadButton(outputId="gtexplot","Download Plot"),
                                                                              actionButton("ex12","Explanation", icon = icon("book-open")),
                                                                              withLoader(verbatimTextOutput("ex_stouffer_gtex"),type="image",loader = "Asset 1.jpg"),

                                                                              br(),br(),

                                                                              h6("One Sided Wilcoxon P-Value"),
                                                                              withLoader(uiOutput("gtex_pstouffer"),type="image",loader = "Asset 1.jpg"),
                                                                              downloadLink('stouffer_table_gtex', 'Download Table', class = "butt",style = "color:yellow"),br(),br(),


                                                                              actionButton("toTop100", "Top", icon = icon("arrow-alt-circle-up")),
                                                               ) )
                                         )


                             )),


                   tabPanel(tags$b("7. Gene Regulatory Network"),
                            tabsetPanel(
                              tabPanel("Cell Cluster",
                                       br(),br(),
                                       sliderInput("grn_obs", "Select Cluster Number ",
                                                   min = 0, max = 20,
                                                   value = 1), br(),
                                       h6("Generate Specific Cluster Table"),br(),


                                       withLoader(verbatimTextOutput("grn_clust"),type="image",loader = "Asset 1.jpg"),

                              ),


                              tabPanel("CellEnrich",

                                       br(),br(),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Transcription Factors and Target Genes)"),
                                      withLoader(simpleNetworkOutput("cell_grn_1"),type="image",loader = "Asset 1.jpg"),
                                      h6("Gene Regulatory Network for CellEnrich (Interaction between Target Genes and Tissue Signatures)"),
                                      withLoader(simpleNetworkOutput("cell_grn"),type="image",loader = "Asset 1.jpg"),


                                       actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                       withLoader(verbatimTextOutput("ex_cellgrn"),type="image",loader = "Asset 1.jpg"),


                                       br(),br(),
                                       h6("Gene Regulatory Table"),
                                       withLoader(uiOutput("cell_grn_table"),type="image",loader = "Asset 1.jpg"),br(),
                                       actionButton("ex_cellgrn","Explanation", icon = icon("book-open")),
                                      downloadLink('download_cell_gt', 'Download Table', class = "butt",style = "color:yellow"),br(),br(),

                                       withLoader(verbatimTextOutput("ex_cellgrn_table"),type="image",loader = "Asset 1.jpg"),
                                       br(), actionButton("toTop6", "Top", icon = icon("arrow-alt-circle-up")),

                              ),tabPanel("TissueEnrich",
                                         tabsetPanel(type = "tabs",
                                                     
                                                     
                                                     
                                                     tabPanel("HPA Data",
                                         br(),br(),

                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("hpa_grn"),type="image",loader = "Asset 1.jpg"),
                                         h6("Gene Regulatory Network for HPA Datasets (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("hpa_grn_1"),type="image",loader = "Asset 1.jpg"),

                                         br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn"),type="image",loader = "Asset 1.jpg"),
                                         br(),br(),

                                         h6("Gene Regulatory Table (HPA)"),
                                         withLoader(uiOutput("hpa_grn_table"),type="image",loader = "Asset 1.jpg"),
                                         br(),
                                         downloadLink('download_hpa_gt', 'Download Table', class = "butt",style = "color:yellow"),br(),br(),
                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_hpagrn_table"),type="image",loader = "Asset 1.jpg"),br(),br(),actionButton("toTop2024", "Top", icon = icon("arrow-alt-circle-up")),

                                         
                                         ),

            
            
            
            tabPanel("GTEx Data",
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Transcription Factors and Target Genes)"),
                                         withLoader(simpleNetworkOutput("gtex_grn_1"),type="image",loader = "Asset 1.jpg"),
                                         h6("Gene Regulatory Network for GTEx data (Interaction between Target Genes and Tissue Signatures)"),
                                         withLoader(simpleNetworkOutput("gtex_grn"),type="image",loader = "Asset 1.jpg"),
                                         br(),


                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),
                                         withLoader(verbatimTextOutput("ex_gtexgrn"),type="image",loader = "Asset 1.jpg"),

                                         br(),br(), h6("Gene Regulatory Table (GTEx)"),
                                         withLoader(uiOutput("gtex_grn_table"),type="image",loader = "Asset 1.jpg"),
                                         br(),

                                         actionButton("ex_grn","Explanation", icon = icon("book-open")),downloadLink('download_gtex_gt', 'Download Table', class = "butt",style = "color:yellow"),
                                         withLoader(verbatimTextOutput("ex_gtexgrn_table"),type="image",loader = "Asset 1.jpg"), br(),br(),
                                         actionButton("toTop7", "Top", icon = icon("arrow-alt-circle-up")),
                             

)
)

 )


                            )

                   )





                    )


             ),










             tabPanel("Tutorial", icon = icon("angle-double-right"),
                     h6( "Check Out the Video Tutorial here",align="center"),
                    tags$iframe(width="1500", height="650", src="EcTracker_Tutorial.mp4", frameborder="0", allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"),
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
