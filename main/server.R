
options(shiny.maxRequestSize = 10000 * 1024 ^ 2)


shinyServer(function(input,output,session){

  session$cache <- cachem::cache_mem(max_size = 100e6)
  session$onSessionEnded(stopApp)

  data_batch<-reactive({
    req(input$file_batch)
    file_b1<-input$file_batch
    inputfile_batch<- read.table(file=file_b1$datapath, header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)

  })


  output$downloadData_batch <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(batch_matrix, file)

    })
output$feature_plot_result<-renderPlot({
  if(input$seurat1==0){
    return()
  }
  FeaturePlot(pb, features = input$feature_plot)
})




  observeEvent(input$testme,{
    if (input$testme) {
      hide('user_input')
    } else if (input$testme)
    {
      hide('user_input')
      hide("testme")
    }
    else {
      show('user_input')
    }
  })
  data <- reactive({
   
                 
    if (input$testme) {
      #write.csv(gse,"gse_new.csv",row.names = FALSE)
      inputfile <- read.csv("gse_new.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)
    }
   # else if (input$testme2) {
   #    # write.csv(gse,"gse_new.csv",row.names = FALSE)
   #    inputfile <- read.csv("pan_exp.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)
   #  }
    
                  
    else if(! is.null(input$File10)){
      req(input$File10)
      File10 <- input$File10
      filesdir = dirname(File10[1,4])

      #inFile = inFile$datapath

      file.rename(File10$datapath[1],paste0(filesdir,'/',File10$name[1]))
      file.rename(File10$datapath[2],paste0(filesdir,'/',File10$name[2]))
      file.rename(File10$datapath[3],paste0(filesdir,'/',File10$name[3]))

      pbmc.data <- Read10X(data.dir = filesdir)
    }

    else {
      # user inputs count matrix
      req(input$file)
      file1 <- input$file
      #if(is.null(file1)){return()}
      inputfile<- read.table(file=file1$datapath, header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)}
                   
  })

data1<-reactive({
  if (input$testme3) {
    #write.csv(meta,"meta.csv",row.names = FALSE)
    inputfile <- read.csv("meta.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)
  }
# else  if (input$testme4) {
#    # write.csv(meta,"meta.csv",row.names = FALSE)
#     inputfile <- read.csv("pan_meta.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)
#   }
 else {
  req(input$meta)

  meta1<-input$meta
inputfile<- read.csv(file=meta1$datapath,sep=",")}

   })

output$merge_meta<-renderPrint({
  if(input$m1==0){
    return()
  }
  cluster_ident<-pb@active.ident

  cluster_ident<-as.data.frame(cluster_ident)
  #names(cluster_ident)[1] <- "Sample"
  write.csv(cluster_ident,file="cluster_ident.csv",row.names = TRUE)

  cluster_ident$Sample<-row.names(cluster_ident)
  cluster_ident
  meta_m<-data1()
 # seurat_m<-as.data.frame(seurat_m)
  #seurat_m<-t(seurat_m)
  #seurat_m<-as.data.frame(seurat_m)
 # seurat_m<-setDT(seurat_m,keep.rownames="...1")
  merger_m<-merge(cluster_ident,meta_m,by.x="Sample",by.y="Sample")
  merger_m
})

meta_file<-reactive({
  cluster_ident<-pb@active.ident

  cluster_ident<-as.data.frame(cluster_ident)
  #names(cluster_ident)[1] <- "Sample"
  write.csv(cluster_ident,file="cluster_ident.csv",row.names = TRUE)

  cluster_ident$Sample<-row.names(cluster_ident)
  cluster_ident
  meta_m<-data1()
  # seurat_m<-as.data.frame(seurat_m)
  #seurat_m<-t(seurat_m)
  #seurat_m<-as.data.frame(seurat_m)
  # seurat_m<-setDT(seurat_m,keep.rownames="...1")
  merger_m<-merge(cluster_ident,meta_m,by.x="Sample",by.y="Sample")
  merger_m
})




output$dirp<-renderPrint({
  if(input$A2==0){
    return()
  }
  req(input$File10)
  File10 <- input$File10
  filesdir = dirname(File10[1,4])

  #inFile = inFile$datapath

  file.rename(File10$datapath[1],paste0(filesdir,'/',File10$name[1]))
  file.rename(File10$datapath[2],paste0(filesdir,'/',File10$name[2]))
  file.rename(File10$datapath[3],paste0(filesdir,'/',File10$name[3]))

  pbmc.data <- Read10X(data.dir = filesdir)
  pb <<- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200,project = "scRNAseq")
  print("Your object has been created, start downstream analysis now!!!!")


})


  output$dat<-renderTable({
    if(is.null(data())){return ()}
    data()[1:10,1:7]
  })

  output$dimen <- renderPrint({
    print("Dimensions")
    dim(data())
  })
  observeEvent(input$seurat1, {
    # Show a modal when the button is pressed

    shinyalert(animation = TRUE,
               imageUrl = "aayu_popup.gif",
               imageWidth = 400,
               imageHeight = 350,confirmButtonText="OK",confirmButtonCol="black")
    #src = "gig.gif"
  })
  observeEvent(input$SCE1, {
    # Show a modal when the button is pressed
    
    shinyalert(animation = TRUE,
               imageUrl = "aayu_popup.gif",
               imageWidth = 400,
               imageHeight = 350,confirmButtonText="OK",confirmButtonCol="black")
    #src = "gig.gif"
  })
  
  observeEvent(input$batch1, {
    # Show a modal when the button is pressed
    
    shinyalert(title ="Please download this file and upload again in Upload matrix box for further analysis.", type = "success")
    #src = "gig.gif"
  })
  
  
  observeEvent(input$action1, {
    # Show a modal when the button is pressed
    
    shinyalert(animation = TRUE,
               imageUrl = "aayu_popup.gif",
               imageWidth = 400,
               imageHeight = 350,confirmButtonText="OK",confirmButtonCol="black")
    #src = "gig.gif"
  })
  observeEvent(input$action11, {
    # Show a modal when the button is pressed
    shinyalert(title = " Markers are being generated... Please check Markers table", type = "success")
  })
output$meta_info<-renderUI({
  if(input$Batch_Correction_method1==0){
    return()
  }
  meta_d<<-data1()
  pb<<-AddMetaData(pb,meta_d)
 choices_meta<<-colnames(pb@meta.data)
 selectInput("meta_header","Select Parameter for Batch Correction",choices = choices_meta)
})






output$meta_info_dropclust<-renderUI({
  if(input$Batch_Correction_method_dropclust1==0){
    return()
  }
  meta_d<-data1()
  colData(sce)<<-DataFrame(meta_d)
  choices_meta_d<<-colnames(sce@colData)
  
  selectInput("meta_header_d","Select parameter for Batch Correction",choices = choices_meta_d)
})




metadata_add<-reactive({
 meta_d<-data1()
  pb<-AddMetaData(pb,meta_d)
})
 
  
  # output$seu_object <- renderPrint({
  #   if(input$seurat1==0){
  #     return()
  #   }
  #   withProgress(message = 'Generating Object...', value = 0,{
  #     pb <<- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200,project = "scRNAseq")
  #     print(pb)
  #     print("Your object has been created, start downstream analysis now!!!!")
  #   })
  # })
  output$seu_object <- renderPrint({
    if(input$seurat1==0){
      return()
    }
    withProgress(message = 'Creating Object...', value = 0,{
      pb <<- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200,project = "scRNAseq")
      print(pb)
      print("Downstream analysis will start now!")
    })
  })
  output$ex_object<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    paste("The object serves as a container that contains both data (like the count matrix) and analyses (like PCA, or clustering results) for a single-cell dataset")
  })

  output$dex_scatter<-renderPrint({
    if(input$SCE1==0){
      return()
    }
    paste("t-SNE representation of the clusters obtained using Louvain algorithm")
  })
  output$dsce_scatter<-renderPrint({
    if(input$SCE1==0){
      return()
    }
    paste("Visualization of cluster assignments for all cells on the t-SNE plot")
  })
  output$dsce_silhouette<-renderPrint({
    if(input$SCE1==0){
      return()
    }
    paste("Check the separateness of the clusters using the silhouette width. Cells with large positive silhouette widths are closer to other cells in the same cluster than to cells in different clusters. Conversely, cells with negative widths are closer to other clusters than to other cells within their own cluster. Each cluster would ideally contain several cells with large positive widths, indicating that it is well-separated from other clusters")
  })
  output$dex_scatter1<-renderPrint({
    if(input$dropclust_batch1==0){
      return()
    }
    paste("Visualization of cells of different batches in t-SNE plot")
  })
  
  output$dex_scatter2<-renderPrint({
    if(input$dropclust_batch1==0){
      return()
    }
    paste("These are per-gene corrected expression values by using their rotation vectors")
  })
  
  
  
  vln<-reactive({
      
    VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

  })
  output$p1<-renderPlot({
    if(input$seurat1==0){
      return()
    }


   VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

  })
  output$ex_variable<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    paste("Visualization of QC metrics as violin plots after filtering the cells based on the number of genes (nFeature_RNA) and reads (nCount_RNA)")
  })

  output$dndPlot <- downloadHandler(
    filename = function(){'VLN_plot.pdf'},
    content = function(file){
      pdf(file)
      print(    vln()
      )
      dev.off()
    })


  observeEvent(input$toTop, {
    js$toTop();
  })
  observeEvent(input$toTop1, {
    js$toTop();
  })
  observeEvent(input$toTop3, {
    js$toTop();
  })
  observeEvent(input$toTop4, {
    js$toTop();
  })
  observeEvent(input$toTop5, {
    js$toTop();
  })
  observeEvent(input$toTop2, {
    js$toTop();
  })
  observeEvent(input$toTop6, {
    js$toTop();
  })
  observeEvent(input$toTop7, {
    js$toTop();
  })
  observeEvent(input$toTop8, {
    js$toTop();
  })
  observeEvent(input$toTop100010, {
    js$toTop();
  })
  
  observeEvent(input$toTop1000101, {
    js$toTop();
  })
  observeEvent(input$toTop10001010, {
    js$toTop();
  })
  observeEvent(input$toTop2021, {
    js$toTop();
  })
  observeEvent(input$toTop2022, {
    js$toTop();
  })
  observeEvent(input$toTop2023, {
    js$toTop();
  })
  observeEvent(input$toTop2024, {
    js$toTop();
  })
  observeEvent(input$toTop2025, {
    js$toTop();
  })
  output$data_norm<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    print("Normalization Done")
  })
  output$ex_normalization<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    paste("It employs a global level normalization method 'LogNormalize' that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor, and log-transforms the result")
  })

  output$variable_features<-renderPrint({
    withProgress(message = 'Creating Feature Plot...', value = 0,{
      
    pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)
      #list_of_variable_features<<-VariableFeatures(pb)
      #list_of_variable_features<<-as.data.frame(list_of_variable_features)
      top10 <<- head(VariableFeatures(pb), 20)
      })
  })


  output$ex_vfeature<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    print("Variable feature helps in subsetting of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). These genes, in downstream analyses, helps to highlight biological signal in single-cell datasets.")
  })


  variableplot<-reactive({
    pb <<- NormalizeData(pb)

    pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)

    top10 <<- head(VariableFeatures(pb), 20)
    plot1 <- VariableFeaturePlot(pb)
    plot2 <- LabelPoints(plot = plot1, points = top10,repel=TRUE,xnudge = 0,ynudge = 0)
    plot(plot2)
  })




  output$plot_vf<-renderPlot({
    if(input$seurat1==0){
      return()
    }
    withProgress(message = 'Creating Variable Feature Plot...', value = 0,{
      
    pb <<- NormalizeData(pb)

    pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)
    top10 <<- head(VariableFeatures(pb), 20)
      plot1 <- VariableFeaturePlot(pb)
      plot2 <- LabelPoints(plot = plot1, points = top10,repel=TRUE,xnudge = 0,ynudge = 0)
     plot(plot2)
    })
  })


  output$vfPlot <- downloadHandler(
    filename = function(){'Vf_plot.pdf'},
    content = function(file){
      pdf(file)
      print(       variableplot()
      )
      dev.off()
    })



  output$scale_data<-renderPrint({
    if(input$seurat1==0){
      return()
    }

    print("Centering and scaling data matrix")
    print("100% Done")
  })
  output$ex_scale<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    print("It is a linear transformation method ('scaling') which is a standard pre-processing step prior to dimensional reduction techniques like PCA. It shifts the expression of each gene so that the mean expression across cells is 0 and scales the expression of each gene, so that the variance across cells is 1. This step gives equal weight to the genes in downstream analyses so that highly-expressed genes do not dominate.")
  })
cca_pca<-reactive({
  
  # if(input$CCA_use){
  #   all.genes <- rownames(pb)
  #   pb <<- ScaleData(pb, features = all.genes)
  #   pb <<- RunPCA(pb, features = VariableFeatures(object = pb),npcs = 50)
  #   pb1 <- subset(pb, cells = colnames(pb)[1:((length(colnames(pb)))/2)])
  #   pb2 <- subset(pb, cells = colnames(x = pb)[(((length(colnames(pb)))/2)+1):length(colnames(pb))])
  #   pb1[["group"]] <- "group1"
  #   pb2[["group"]] <- "group2"
  #   pb<<- RunCCA(object1 = pb1, object2 = pb2)
  #   pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)
  #   pb <<- ScaleData(pb, features = all.genes)
  #   pb <<- RunPCA(pb, features = VariableFeatures(object = pb),npcs = 50)
  #   RunPCA(pb, reduction.use = "cca.aligned")
  #   DimPlot(pb, reduction = "cca", pt.size = 4, group.by = "ident")
  # }
  # else if(input$SCT_use)
  # {
  #   pb<<- PercentageFeatureSet(pb, pattern = "^MT-", col.name = "percent.mt")
  #   pb <<- SCTransform(pb, vars.to.regress = "percent.mt", verbose = FALSE)
  #   pb <<- RunPCA(pb, verbose = FALSE)
  #   pb <<- RunUMAP(pb, dims = 1:30, verbose = FALSE)
  #   
  #   pb <<- FindNeighbors(pb, dims = 1:30, verbose = FALSE)
  #   pb <<- FindClusters(pb, verbose = FALSE)
  #   DimPlot(pb, label = TRUE)
  # }
  # else {
    #  withProgress(message = 'Scaling and PCA ...', value = 0,{
    all.genes <- rownames(pb)
    pb <<- ScaleData(pb, features = all.genes)
    pb <<- RunPCA(pb, features = VariableFeatures(object = pb),npcs = 50)
    DimPlot(pb, reduction = "pca", pt.size = 4, group.by = "ident")
    # }
  #})
})


    output$pca<-renderPlot({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Generating PCA Plot...', value = 0,
                   {
    cca_pca()
                   })
  })

    # pcaplot<-reactive({
    #   
    #   all.genes <- rownames(pb)
    #   pb <<- ScaleData(pb, features = all.genes)
    #   pb <<- RunPCA(pb, features = VariableFeatures(object = pb),npcs = 48)
    #   DimPlot(pb, reduction = "pca", pt.size = 4, group.by = "ident")
    # })
    
    
  output$pcaPlot <- downloadHandler(
      filename = function(){'PCA_plot.pdf'},
      content = function(file){
        pdf(file)
        print(        cca_pca()

        )
        dev.off()
      })

  output$ex_pca<-renderPrint({
    if(input$seurat1==0){
      return()
    }
    paste("It performs linear dimensionality reduction")
  })
  output$hex_pca<-renderPrint({
    if(input$Batch_Correction_method==0){
      return()
    }
    paste("Plot shows how well the datasets are integrated after Harmony")
  })
  output$hex_harmplot<-renderPrint({
    if(input$Batch_Correction_method==0){
      return()
    }
    paste("Seurat object, updated with the corrected Harmony coordinates")
  })
  
   heatmapplot<-reactive({
     
    DimHeatmap(pb, dims = 1:6, cells = 500, balanced = TRUE)

 })


    output$heatmap<-renderPlot({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Generating Heatmap...', value = 0,{
        
      DimHeatmap(pb, dims = 1:6, cells = 500, balanced = TRUE)
        })

    })
    output$ex_heatmap<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      paste("It allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores.")
    })
    output$ex_karyo<-renderPrint({
      if(input$seurat1==0){
        return()
      }
 paste("It gives the chromosomal location of the selected gene.")
    })
    output$ex_fold_change<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Bar plot represents the fold-change values of the tissue-specific gene enrichment")
    })
   
    
    output$ex_log_p_value<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Bar plot represents the tissue-specific gene enrichment (???Log10(P???Value))")
    })
    
     output$ex_feature<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      paste("Visualization of single cells on dimensional reduction plot based on a user selected gene")
    })
    
    
    output$heatmapPlot <- downloadHandler(
      filename = function(){'HEATMAP_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          heatmapplot())
        dev.off()
      })

    jack<-reactive({  
      pb <<- JackStraw(pb, num.replicate = 100)
      pb <<- ScoreJackStraw(pb, dims = 1:20)
      JackStrawPlot(pb, dims = 1:15)
    
      
    })
    output$jackstraw<-renderPlot({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Generating JackStraw Plot...', value = 0,{
        
      pb <<- JackStraw(pb, num.replicate = 100)
      pb <<- ScoreJackStraw(pb, dims = 1:20)
      #pdf(file="Main_pipeline/GSE756881/Jackstraw_plot.pdf")
      JackStrawPlot(pb, dims = 1:15)
      })
      
    })
    output$ex_jackstraw<-renderPrint({
      if(input$seurat1==0){
        return()
      }
       paste("The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 'Significant' PCs will show strong enrichment of features with low p-values (solid curve above the dashed line)")
          })

    output$jackstrawPlot <- downloadHandler(
      filename = function(){'Jackstraw_plot.pdf'},
      content = function(file){
        pdf(file)
        print(jack()

          )
        dev.off()
      })

    output$ex_elbow<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      paste("A ranking of principal components based on the percentage of variance explained by each of them")
    })
    output$elbow<-renderPlot({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Generating Elbow Plot...', value = 0,{
        
      ElbowPlot(pb)
})
    })


    elbow<-reactive({
      ElbowPlot(pb)
    })
    output$elbowPlot <- downloadHandler(
      filename = function(){'elbow_plot.pdf'},
      content = function(file){
        pdf(file)
        print(       elbow()

        )
        dev.off()
      })


    output$neigh<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Clustering...', value = 0,{
        
      pb <<- FindNeighbors(pb, dims = 1:10)
      pb <<- FindClusters(pb, resolution = 0.5)

       head(Idents(pb), 5)
})
    })
    output$ex_umap<-renderPrint({
      if(input$seurat1==0){
        return()
      }
paste("The goal of this algorithm is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimensionality reduction plots")
          })
    output$hex_umap<-renderPrint({
      if(input$Batch_Correction_method==0){
        return()
      }
      paste("The goal of this algorithm is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimensionality reduction plots")
    })
    
    output$ex_umapp<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      paste("The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots.")

      })
    output$ex_stouffer<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Box plot depicting Stouffer's score computed based on the interaction between related genes across indicated cell types.")
    })


    output$ex_stouffer_hpa<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Box plot depicting Stouffer's score computed based on the interaction between related genes across indicated cell types.")
    })
    output$ex_stouffer_gtex<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Box plot depicting Stouffer's score computed based on the interaction between related genes across indicated cell types.")
    })

    output$ex_hpagrn<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Gene regulatory method helps users to understand which transcription factor helps the target gene in expressing a particular signature. Two graphs are plotted which shows the interaction between transcription factors and target genes (Green links) and the other between target genes and signatures (Pink links)")
    })
    output$ex_hpagrn_table<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("The table represents the target genes with their respective transcription factors and signatures. Here, mor stands for mode of regulation where -1, +1 are the two modes, while confidence level is represented by A, B and C where A means higher confidence.") 
    })
    output$ex_gtexgrn<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Gene regulatory method helps users to understand which transcription factor helps the target gene in expressing a particular signature. Two graphs are plotted which shows the interaction between transcription factors and target genes (Green links) and the other between target genes and signatures (Pink links)")
    })
    output$ex_gtexgrn_table<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("The table represents the target genes with their respective transcription factors and signatures. Here, mor stands for mode of regulation where -1, +1 are the two modes, while confidence level is represented by A, B and C where A means higher confidence.")
    })
    output$ex_cellgrn<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("Gene regulatory method helps users to understand which transcription factor helps the target gene in expressing a particular signature. Two graphs are plotted which shows the interaction between transcription factors and target genes (Green links) and the other between target genes and signatures (Pink links)")
    })
    output$ex_cellgrn_table<-renderPrint({
      if(input$action1==0){
        return()
      }
      paste("The table represents the target genes with their respective transcription factors and signatures. Here, mor stands for mode of regulation where -1, +1 are the two modes, while confidence level is represented by A, B and C where A means higher confidence.")
    })


    marker_display_option<-reactive({
      if (input$test_use=="Wilcoxon")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "wilcox")}
      
      else if (input$test_use=="Bimod (Likelihood-ratio test for single cell gene expression)")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "bimod")}
      else if (input$test_use=="ROC analysis")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "roc")}
      
      else if (input$test_use=="Student's t-test")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "t")}
      
      else if (input$test_use=="Negative binomial generalized linear model")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "negbinom")}
      else if (input$test_use=="Poisson generalized linear model")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "poisson")}
      
      else if (input$test_use=="Logistic Regression")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "LR")}
      
      else if (input$test_use=="MAST")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "MAST")}
      else if (input$test_use=="DESeq2")
      { pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc,test.use = "DESeq2")}
      
    })
    
    output$markers_display<-renderTable({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Marker Identification...', value = 0,{
     marker_display_option()
   
      # identity_cell<<-Idents(pb)

       write.table(pb.markers,file="Markers_info.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)
      head(pb.markers)
     })
    })


     #coding for the tissue enrich

    tissue<-reactive({
      clus_assay<-GetAssayData(pb)

      clus_assay<-as.data.frame(clus_assay)

      # write.csv(clus_assay,file="clus_assay.csv",row.names = TRUE)
    })

    cluster_ident_ <-reactive({
      cluster_ident<-pb@active.ident

      cluster_ident<-as.data.frame(cluster_ident)
      # names(cluster_ident)[1] <- "Sample"
      write.csv(cluster_ident,file="cluster_ident.csv",row.names = TRUE)
      cluster_ident
    })

    output$cluster_ident<- renderPrint(
      {
        if(input$seurat1==0){
          return()
        }
        tissue()####### after added
        cluster_ident_()


      })
    output$downloadData_cluster_ident <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(cluster_ident_(), file)

      })


    output$downloadData_meta <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(meta_file(), file)

      })

    marker_table<- reactive({

      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      #  write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      # marker_dataframe

    })
    marker_table_grn<- reactive({
      ################# adding for tissue enrich

      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

      gene_name_cluster<- marker_dataframe$gene
      # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe

    })

    output$cluster_name<-renderPrint({
      if(input$action24==0){
        return()
      }
      marker_table()


    })

    output$grn_clust<-renderPrint({
      if(input$action1==0){
        return()
      }
      marker_table_grn()
    })

    output$downloadData_cluster_name <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(marker_table(), file)

      })

    tissue_enrich<-reactive({


      #gene_name_cluster<- marker_dataframe$gene

      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      #write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe

      inputGenes<-gene_name_cluster
      inputGenes

    })

    output$tissue_detail<-renderPrint({
      if(input$action1==0){
        return()
      }
      #marker_table()
      withProgress(message = 'Starting Tissue Enrichment...', value = 0,
                   {
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)



      tab<- enrichmentOutput[,-4]
      tab
      log_df<-subset(tab, tab$Log10PValue!=0.00000000 & tab$fold.change!=0.00000000)
      log_df
})
    })


    output$downloadtissuedata <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(enrichmentOutput, file)

      })

    plot_fold_t<-reactive({
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      # marker_dataframe
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)


      tab<- enrichmentOutput[,-4]
      tab
      log_df<-subset(tab, tab$Log10PValue!=0.00000000 & tab$fold.change!=0.00000000)
      log_df

      ggplot(log_df,aes(x=reorder(Tissue,-fold.change),y=fold.change,label = Tissue.Specific.Genes,fill = Tissue))+
        geom_bar(stat = 'identity', width=0.5)+
        labs(x='', y = 'Fold change')+
        theme_bw()+
        theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

    })

    output$plot_fold<-renderPlot({
      if(input$action1==0){
        return()
      }
      withProgress(message = 'Calculating Fold Change...', value = 0,
                   {
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      tab<- enrichmentOutput[,-4]
      tab
      log_df<-subset(tab, tab$Log10PValue!=0.00000000 & tab$fold.change!=0.00000000)
      log_df

      ggplot(log_df,aes(x=reorder(Tissue,-fold.change),y=fold.change,label = Tissue.Specific.Genes,fill = Tissue))+
        geom_bar(stat = 'identity', width=0.5)+
        labs(x='', y = 'Fold change')+
        theme_bw()+
        theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
})

    })

    output$tissue_fold1 <- downloadHandler(
      filename = function(){'fold_change.pdf'},
      content = function(file){
        pdf(file)
        print(
          plot_fold_t()

        )
        dev.off()
      })


    output$tissue_log1 <- downloadHandler(
      filename = function(){'log_pvalue.pdf'},
      content = function(file){
        pdf(file)
        print(

          logfold()


        )
        dev.off()
      })

    output$plot_log<-renderPlot({
      if(input$action1==0){
        return()
      }
      withProgress(message = 'Calculating Log P-Value...', value = 0,
                   {
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe

      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      tab<- enrichmentOutput[,-4]
      tab
      log_df<-subset(tab, tab$Log10PValue!=0.00000000 & tab$fold.change!=0.00000000)
      log_df
      ggplot(log_df,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
        geom_bar(stat = 'identity', width=0.5)+
        labs(x='', y = '-LOG10(P-Adjusted)')+
        theme_bw()+
        theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
                   })
        })

    logfold<-reactive({

      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe


      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)


      tab<- enrichmentOutput[,-4]
      tab
      log_df<-subset(tab, tab$Log10PValue!=0.00000000 & tab$fold.change!=0.00000000)
      log_df
      ggplot(log_df,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
        geom_bar(stat = 'identity', width=0.5)+
        labs(x='', y = '-LOG10(P-Adjusted)')+
        theme_bw()+
        theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())})

    output$kay_name<-renderPrint({


      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      marker_dataframe
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)



      enrichmentOutput

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])

      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      seGroupInf<-output2[[3]][[input$kar1]]
      groupInf<-data.frame(assay(seGroupInf))

      groupInf$Gene


    })


    umap_plot<-reactive({
      DimPlot(pb, reduction = "umap",pt.size = 3)

    })

    output$umap<-renderPlot({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Generating UMAP...', value = 0,{
        
      pb <<- FindNeighbors(pb, dims = 1:10)
      pb <<- FindClusters(pb, resolution = 0.5)

      head(Idents(pb), 5)
      pb <<- RunUMAP(pb, dims = 1:10)
      umap_mat<<-pb[["umap"]]@cell.embeddings
      umap_mat<<-as.data.frame(umap_mat)
      DimPlot(pb, reduction = "umap",pt.size = 3)
})
    })
    output$umapPlot <- downloadHandler(
      filename = function(){'UMAP_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          umap_plot()


        )
        dev.off()
      })

    list<-reactive({
      lst1=list()
      for (i in a){lst[[i]]=assay(seGroupInf[[i]])}
      lst1
    })
    pbmarker<-reactive({
      pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)

    })

    output$markers<-renderPrint({
      if(input$seurat1==0){
        return()
      }
      withProgress(message = 'Marker Identification...', value = 0,{
        
      pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)
       write.table(pb.markers,file="Markers_info.csv",
                 sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)
       print("Markers identified!!!!!")
       })
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(pb.markers, file)

      })



    output$stouffer_sc_plot<-renderPlot({
      if(input$action1==0){
        return()
      }
      withProgress(message = 'Generating Stouffer Plot...', value = 0,
                   {
        matric<-pb@assays$RNA@counts

        cells_sum <- Matrix::colSums(matric>=3)

        matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

        matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

        cells_sum<-Matrix::rowSums(t(matric))

        matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
        r_name<-rownames(matric)

       # Stouffer_CellEnrich()
        media_genes<-abc[[input$st_sc]]@geneIds
        media_genes
        same_genes<-intersect(r_name,media_genes)
        if(length(same_genes)==1)
          message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
        if(length(same_genes)==0)
          message="no common gene "
        if(length(same_genes)>1){
          matric= matric[same_genes,]####### didnot work why?

          #Further insuring have cells with atleast 10% expressed genes
          matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
          print(dim(matric))
          seurat_RNA_matric<-matric
          print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
          print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
          print(dim(seurat_RNA_matric))
          ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
          seurat_RNA_matric[seurat_RNA_matric==0]=1
          seurat_RNA_matric<-log2(seurat_RNA_matric)
          seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
          ###Stouffer score and then boxplot for each cell type
          stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
          print(sum(is.na(stouffer_score)))
          cell_types<-Idents(pb)[names(stouffer_score)]
          unique_cell_types=unique(cell_types)
          Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
          ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
            geom_boxplot(position=position_dodge(.2)) +
            geom_jitter(shape=16, position=position_jitter(.1)) +
            theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}

        })
    })













    output$de_stouffer<-renderPlot({
      if(input$action230==0){
        return()
      }
      stuouffer()
    })

    output$stoufferPlot <- downloadHandler(
      filename = function(){'Stouffer_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
            geom_boxplot(position=position_dodge(.2)) +
            geom_jitter(shape=16, position=position_jitter(.1)) +
            theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            scale_fill_brewer(palette="Dark2")

        )
        dev.off()
      })



    output$de_karyotype<-renderPlot({
      if(input$seurat1==0){
        return()
      }

      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      dm.genes <- genes(txdb)
      dm.genes

      library("org.Hs.eg.db")
      kar<-pb.markers
      gene<-(input$feature_plot)
      gene
      gene_id<-unlist(mget(x=gene,envir=org.Hs.egALIAS2EG))
      gene_id

      #Merging the data frames
      gened<-as.data.frame(gene_id)
      #gened$name<-row.names(gened)
      askk<-data.frame("gene"=rownames(gened))
      askk$p_val<-kar[gene,]$p_val
      askk$avg_log2FC<-kar[gene,]$avg_log2FC
      askk$pct.1<-kar[gene,]$pct.1
      askk$pct.2<-kar[gene,]$pct.2
      askk$p_val_adj<-kar[gene,]$p_val_adj
      askk$cluster<-kar[gene,]$cluster
      askk$gene_id<-gened$gene_id
      rownames(askk)<-askk$gene_id



      #ask<-merge(kar,gened,by.x="gene",by.y="name")

      ##Final data frame with start and end
      mcols(dm.genes) <- askk[names(dm.genes), c("avg_log2FC","p_val","p_val_adj","cluster","gene")]
      head(dm.genes, n=4)

      ##
      library(karyoploteR)
      filtered.dm.genes_mat <- dm.genes[!is.na(dm.genes$p_val_adj)]
      filtered.dm.genes<-filtered.dm.genes_mat[gened$gene_id,]
      ord<-as.data.frame(filtered.dm.genes)
      rownames(ord)<-ord$gene
      log.pval <- -log10(filtered.dm.genes$p_val_adj)
      mcols(filtered.dm.genes)$log.pval <- log.pval
      filtered.dm.genes

      range(filtered.dm.genes$avg_log2FC)
      fc.ymax <- ceiling(max(abs(range(filtered.dm.genes$avg_logFC))))
      fc.ymin <- -fc.ymax

      cex.val <- sqrt(filtered.dm.genes$log.pval)/2

      col.over <- "#00A6EDAA"
      sign.col <- rep(col.over, length(filtered.dm.genes))
      chrs<-unique(ord$seqnames)

      ##Final Ploting
      kp <- plotKaryotype(genome="hg19")
      #kpPoints(kp, data=filtered.dm.genes, y=filtered.dm.genes$avg_log2FC,  cex=cex.val, ymax=fc.ymax, ymin=fc.ymin,r1=0.1,col=sign.col)
      #kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin,r1=0.8)
      #kpPlotMarkers(kp, filtered.dm.genes, labels = rownames(ord), text.orientation = "horizontal", ignore.chromosome.ends=TRUE, r0=0.01,marker.parts=c(0.1,0.1,0.1),label.dist=0.03)
      #kpAddLabels(kp, labels = "log2FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin,r1=0.8)
      #gene.mean <- start(filtered.dm.genes) + (end(filtered.dm.genes) - start(filtered.dm.genes))/2
      #kpSegments(kp, chr=as.character(seqnames(filtered.dm.genes)), x0=gene.mean, x1=gene.mean, y0=filtered.dm.genes$avg_log2FC, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=0.8)

      ##Selected Genes
      gen<-filtered.dm.genes$gene

      # gene_list<- c(input$detected_tissue)
      col.under <- "#0000FF"
      sign.col[filtered.dm.genes$gene==gen] <- col.under
      kpPoints(kp, data=filtered.dm.genes, y=filtered.dm.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin,r1=0.1,col=sign.col)
      #selected.dm.genes <- filtered.dm.genes[filtered.dm.genes$gene]
      kpPlotMarkers(kp, filtered.dm.genes, labels = gen, text.orientation = "horizontal", ignore.chromosome.ends=TRUE, r0=0.01,marker.parts=c(0.1,0.1,0.1),label.dist=0.005,label.color =  "#0000FF")







    })
    output$karyoplott <- downloadHandler(
      filename = function(){'karyoplot_plot.pdf'},
      content = function(file){
        pdf(file)
        print(    kayro())

        dev.off()
      })

  kayro<-reactive({
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    dm.genes <- genes(txdb)
    dm.genes

    library("org.Hs.eg.db")
    kar<-pb.markers
    gene<-(input$feature_plot)
    gene
    gene_id<-unlist(mget(x=gene,envir=org.Hs.egALIAS2EG))
    gene_id

    #Merging the data frames
    gened<-as.data.frame(gene_id)
    #gened$name<-row.names(gened)
    askk<-data.frame("gene"=rownames(gened))
    askk$p_val<-kar[gene,]$p_val
    askk$avg_log2FC<-kar[gene,]$avg_log2FC
    askk$pct.1<-kar[gene,]$pct.1
    askk$pct.2<-kar[gene,]$pct.2
    askk$p_val_adj<-kar[gene,]$p_val_adj
    askk$cluster<-kar[gene,]$cluster
    askk$gene_id<-gened$gene_id
    rownames(askk)<-askk$gene_id



    #ask<-merge(kar,gened,by.x="gene",by.y="name")

    ##Final data frame with start and end
    mcols(dm.genes) <- askk[names(dm.genes), c("avg_log2FC","p_val","p_val_adj","cluster","gene")]
    head(dm.genes, n=4)

    ##
    library(karyoploteR)
    filtered.dm.genes_mat <- dm.genes[!is.na(dm.genes$p_val_adj)]
    filtered.dm.genes<-filtered.dm.genes_mat[gened$gene_id,]
    ord<-as.data.frame(filtered.dm.genes)
    rownames(ord)<-ord$gene
    log.pval <- -log10(filtered.dm.genes$p_val_adj)
    mcols(filtered.dm.genes)$log.pval <- log.pval
    filtered.dm.genes

    range(filtered.dm.genes$avg_log2FC)
    fc.ymax <- ceiling(max(abs(range(filtered.dm.genes$avg_logFC))))
    fc.ymin <- -fc.ymax

    cex.val <- sqrt(filtered.dm.genes$log.pval)/2

    col.over <- "#00A6EDAA"
    sign.col <- rep(col.over, length(filtered.dm.genes))
    chrs<-unique(ord$seqnames)

    ##Final Ploting
    kp <- plotKaryotype(genome="hg19")
    #kpPoints(kp, data=filtered.dm.genes, y=filtered.dm.genes$avg_log2FC,  cex=cex.val, ymax=fc.ymax, ymin=fc.ymin,r1=0.1,col=sign.col)
    #kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin,r1=0.8)
    #kpPlotMarkers(kp, filtered.dm.genes, labels = rownames(ord), text.orientation = "horizontal", ignore.chromosome.ends=TRUE, r0=0.01,marker.parts=c(0.1,0.1,0.1),label.dist=0.03)
    #kpAddLabels(kp, labels = "log2FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin,r1=0.8)
    #gene.mean <- start(filtered.dm.genes) + (end(filtered.dm.genes) - start(filtered.dm.genes))/2
    #kpSegments(kp, chr=as.character(seqnames(filtered.dm.genes)), x0=gene.mean, x1=gene.mean, y0=filtered.dm.genes$avg_log2FC, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=0.8)

    ##Selected Genes
    gen<-filtered.dm.genes$gene

    # gene_list<- c(input$detected_tissue)
    col.under <- "#0000FF"
    sign.col[filtered.dm.genes$gene==gen] <- col.under
    kpPoints(kp, data=filtered.dm.genes, y=filtered.dm.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin,r1=0.1,col=sign.col)
    #selected.dm.genes <- filtered.dm.genes[filtered.dm.genes$gene]
    kpPlotMarkers(kp, filtered.dm.genes, labels = gen, text.orientation = "horizontal", ignore.chromosome.ends=TRUE, r0=0.01,marker.parts=c(0.1,0.1,0.1),label.dist=0.005,label.color =  "#0000FF")




  })






    output$f_or<-renderUI({
      or<-read.csv("final_list.csv")
      row_names<-row.names(pb)
      row_names<-as.data.frame(row_names)
      colnames(row_names)<-"Receptor_names"
      factor0<<-merge(row_names,or,by.x="Receptor_names",by.y="Symbol")
      factor0
      selectInput("detected_or","Select Receptor",choices = factor0)
      #paste("Number of ORs detected",nrow(factor0))
    })

    output$or_umap<-renderPlot({
      FeaturePlot(pb, features = c(input$detected_or))
    })


    output$orPlot <- downloadHandler(
      filename = function(){'OR_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          FeaturePlot(pb, features = c(input$detected_or))


        )
        dev.off()
      })


    output$t_or<-renderUI({
      if(input$action1==0){
        return()
      }
      withProgress(message = 'Finding Tisssue Specific Genes...', value = 0,{
                   
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

      gene_name_cluster<- marker_dataframe$gene
      # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)


      #tissue_subset <-subset(enrichmentOutput$Tissue, enrichmentOutput$fold.change!= 0)


      # tissue_subset <-as.data.frame(tissue_subset)

      enrichmentOutput

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])

      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      zzz<-  enrichmentOutput$Tissue
      zzz


      selectInput("name_tissue","Select Tissue",choices = c(zzz))
})

    })

    gene_tissue_detail<-reactive({
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

    gene_name_cluster<- marker_dataframe$gene
    # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
    gene_name_cluster
      #Retrieval of input tissue-specific genes
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)


      #tissue_subset <-subset(enrichmentOutput$Tissue, enrichmentOutput$fold.change!= 0)


      # tissue_subset <-as.data.frame(tissue_subset)

      enrichmentOutput

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])

      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      seGroupInf<-output2[[3]][[input$name_tissue]]
      groupInf<-data.frame(assay(seGroupInf))
      write.csv(groupInf,"tissue_sp_gene.csv")

    })
   output$tt_or<-renderUI({
     if(input$action1==0){
       return()
     }
    # gene_tissue_detail()
     marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$obs )

     gene_name_cluster<- marker_dataframe$gene
     # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
     gene_name_cluster

    inputGenes <-gene_name_cluster
     inputGenes

     gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
     output2<-teEnrichment(gs, rnaSeqDataset = input$dataset)#here
     output2
     enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

     write.table(enrichmentOutput,file="enrichment_output.csv",
                 sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

     seEnrichmentOutput<-output2[[1]]

     enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


     enrichmentOutput$Tissue<-row.names(enrichmentOutput)


     #tissue_subset <-subset(enrichmentOutput$Tissue, enrichmentOutput$fold.change!= 0)


     # tissue_subset <-as.data.frame(tissue_subset)

     enrichmentOutput

     seEnrichmentOutput<-output2[[1]]

     enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])

     enrichmentOutput$Tissue<-row.names(enrichmentOutput)
     seGroupInf<-output2[[3]][[input$name_tissue]]
     groupInf<-data.frame(assay(seGroupInf))
     write.csv(groupInf,"tissue_sp_gene.csv")
     ttor<- read.csv("tissue_sp_gene.csv")
     # ttor<-as.data.frame(ttor)
     #ttor<-ttor$Gene
     row_names1<-row.names(pb)
     row_names1<-as.data.frame(row_names1)
    colnames(row_names1)<-"Receptor_names"
     factor100<<-merge(row_names1,ttor,by.x="Receptor_names",by.y="Gene")
     factor1000<-factor100$Receptor_names
     factor1000

     selectInput("detected_tissue","Select Transcript",choices = factor1000)
     #paste("Number of ORs detected",nrow(factor0))
   })

    # output$tor_umap<-renderPlot({
    #   FeaturePlot(pb, features = c(input$detected_tissue))
    # })


    output$torPlot <- downloadHandler(
      filename = function(){'umap_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          FeaturePlot(pb, features = c(input$feature_plot))


        )
        dev.off()
      })



    output$user_umap<-renderPlot({
      FeaturePlot(pb, features = c(input$usergene))
    })
    output$uorPlot <- downloadHandler(
      filename = function(){'UOR_plot.pdf'},
      content = function(file){
        pdf(file)
        print(
          FeaturePlot(pb, features = c(input$usergene))


        )
        dev.off()
      })




    output$grn<- renderTable({
      if(input$action6000==0){
        return()
      }


      #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc )### have to make chnage in threshold

      ##########************************* read dorothea_regulon file
      dorothea_regulon_human <-doro

      ## We obtain the regulons based on interactions with confidence level A, B and C
      regulon <- dorothea_regulon_human %>%
        dplyr::filter(confidence %in% c("A","B","C"))

      ## We compute Viper Scores
      pb <- run_viper(pb, regulon,
      )

      ## We compute the Nearest Neighbours to perform cluster
      DefaultAssay(object = pb) <- "dorothea"
      pb <- ScaleData(pb)
      pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
      pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
      pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

      pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

      pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = input$log_fc, verbose = FALSE)

      DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

      ## We transform Viper scores, scaled by seurat, into a data frame to better

      ## handling the results
      viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                      assay = "dorothea") %>%
        data.frame() %>%
        t()

      ## We create a data frame containing the cells and their clusters
      CellsClusters <- data.frame(cell = names(Idents(pb)),
                                  cell_type = as.character(Idents(pb)),
                                  stringsAsFactors = FALSE)


      ################## my addition ###################################
      viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
      viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
      viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

      summarized_viper_scores <- viper_scores_clusters %>%
        group_by(variable, cell_type) %>%
        summarise(avg = mean(value),
                  std = sd(value))
      ## We select the 20 most variable TFs. (20*9 populations = 180)

      highly_variable_tfs <- summarized_viper_scores %>%
        group_by(variable) %>%
        mutate(var = var(avg))  %>%
        ungroup() %>%
        top_n(180, var) %>%
        distinct(variable)

      ## We prepare the data for the plot
      summarized_viper_scores_df <- summarized_viper_scores %>%
        semi_join(highly_variable_tfs, by = "variable") %>%
        dplyr::select(-std) %>%
        spread(variable, avg) %>%
        data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

      transpose_summarise_df<-t(summarized_viper_scores_df)

      merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
      final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
      marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

      gene_name_cluster<- marker_dataframe$gene
      # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
      gene_name_cluster
      inputGenes <-gene_name_cluster
      inputGenes

      gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())# changes TE custom to te
      output2<-teEnrichment(gs)#here
      output2
      enrichmentOutput<<-setNames(data.frame(assay(output2[[1]]), row.names = rowData(output2[[1]])[,1]), colData(output2[[1]])[,1])

      write.table(enrichmentOutput,file="enrichment_output.csv",
                  sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])


      enrichmentOutput$Tissue<-row.names(enrichmentOutput)


      #tissue_subset <-subset(enrichmentOutput$Tissue, enrichmentOutput$fold.change!= 0)


      # tissue_subset <-as.data.frame(tissue_subset)

      enrichmentOutput

      seEnrichmentOutput<-output2[[1]]

      enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])

      enrichmentOutput$Tissue<-row.names(enrichmentOutput)

      seGroupInf<-output2[[3]]
      m=matrix(0,nrow=35,ncol= 20)
      rownames(m)=names(seGroupInf)
      a=names(seGroupInf)

      for (i in a){if(length(assay(seGroupInf[[i]])[,"Gene"])!=0)m[i,1:length(assay(seGroupInf[[i]])[,"Gene"])]=assay(seGroupInf[[i]])[,"Gene"]}
      m<-as.data.frame(m)
      names(m)[1:10]<- "target"
      tm<-t(m)
      melted_m<- melt(tm,id="tissuenames",measure.vars = c("target","target.1","target.2","target.3","target.4","target.5","target.6","target.7","target.8","target.9","target.10"))
      melted_m<-unique(melted_m)
      updated_myData <- subset(melted_m, value!= 0)
      finaldata_tissuenerich<-merge(final_data_plot,updated_myData,by.x="target", by.y="value")
      colnames(finaldata_tissuenerich)[which(names(finaldata_tissuenerich) == "Row.names")] <- "TF"

      final_data_toshow<- finaldata_tissuenerich  %>% dplyr::select(target,Var2,TF,mor)
      colnames(final_data_toshow)[which(names(final_data_toshow) == "Var2")] <- "Tissue"
      colnames(final_data_toshow)[which(names(final_data_toshow) == "target")] <- "Target"

      final_data_toshow
      #chordDiagram(finaldata_tissuenerich)
      #final_data_toshow$mor <- as.factor(final_data_toshow$mor)

     # ggplot(as.data.frame(final_data_toshow),
     #        aes(axis1 = Tissue, axis2 = TF, axis3 = Target)) +
     ###   geom_alluvium(aes(fill = mor), width = 1/12) +
       ## geom_stratum(width = 1/12, fill = "black", color = "grey") +
       # scale_x_continuous(breaks = 1:3, labels = c("Tissue", "TF", "Target")) +
        #scale_fill_brewer(type = "qual", palette = "Set1") +
        #ggtitle("Gene regulatory network") +
        #coord_flip() +
        #ggfittext::geom_fit_text(aes(label = after_stat(stratum)), reflow = TRUE, stat = "stratum", width = 1/4, min.size = 1, colour = "white", size = 10) +
        #theme_minimal()



    })
      # alluvial::alluvial(finaldata_tissuenerich,freq = dim(finaldata_tissuenerich),col = "blue")

      #chordDiagram(finaldata_tissuenerich)





    output$grnPlot <- downloadHandler(
      filename = function(){'GRN.pdf'},
      content = function(file){
        pdf(file)
        print(   plot(plotg, layout=-l$layout[,2:1],vertex.size = 25, vertex.color = "PINK",
                      vertex.frame.color = "RED", vertex.label.color = "black",
                      vertex.label.cex = 1)


        )
        dev.off()
      })


    AUCell_table<-reactive({
      au_mat<-pb@assays$RNA@counts
      au_mat<-as.matrix(au_mat)

      cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
      cells_rankings
      #save(cells_rankings, file="cells_rankings.RData")
      cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
      #save(cells_AUC, file="cells_AUC.RData")
      set.seed(123)
      par(mfrow=c(1,1))
      cells_assignment <- AUCell_exploreThresholds(cells_AUC, assign=TRUE)
      score_cell<-as.matrix(cells_assignment)
   score_cell

    })
    AUCell<-reactive({
      au_mat<-pb@assays$RNA@counts
      au_mat<-as.matrix(au_mat)

      cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
      cells_rankings
      #save(cells_rankings, file="cells_rankings.RData")
      cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
      #save(cells_AUC, file="cells_AUC.RData")
      set.seed(123)
      par(mfrow=c(1,1))
      cells_assignment <- AUCell_exploreThresholds(cells_AUC, assign=TRUE)
      df2 <- as.data.frame(lapply(cells_assignment, unlist))
      sumByGene <- apply(au_mat, 1, sum)
      exprMatSubset <- au_mat[which(sumByGene>0),]
      logMatrix <- log2(exprMatSubset+1)

     # umap_mat<-pb[["umap"]]@cell.embeddings
      umap_mat<-as.data.frame(umap_mat)
      selectedThresholds <- getThresholdSelected(cells_assignment)

      for(geneSetName in names(selectedThresholds))
      {
        nBreaks <- 5 # Number of levels in the color palettes
        # Color palette for the cells that do not pass the threshold
        colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
        # Color palette for the cells that pass the threshold
        colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)

        # Split cells according to their AUC value for the gene set
        passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
        if(sum(passThreshold) >0 )
        {
          aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)

          # Assign cell color
          cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
                         setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))

          # Plot
          plot(umap_mat, main=geneSetName,
               sub="Pink/red cells pass the threshold",
               col=cellColor[rownames(umap_mat)], pch=16)
        }
      }



    })

#       output$au_cell<-renderPlot({
#     if(input$action1==0){
#       return()
#     }
#         withProgress(message = 'Cell Assignments', value = 0,
#                      {
#         au_mat<-pb@assays$RNA@counts
#         au_mat<-as.matrix(au_mat)
# 
#         cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
#         cells_rankings
#        # save(cells_rankings, file="cells_rankings.RData")
#         cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
#        # save(cells_AUC, file="cells_AUC.RData")
#         set.seed(123)
#         par(mfrow=c(1,1))
#         cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
# })
#         
#          })

      output$au_cell_R<-renderPlot({
        if(input$action1==0){
          return()
        }
        withProgress(message = 'Cell Assignments...', value = 0,
                     {
                       au_mat<-pb@assays$RNA@counts
                       au_mat<-as.matrix(au_mat)
                       
                       cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
                       cells_rankings
                       # save(cells_rankings, file="cells_rankings.RData")
                       cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
                       # save(cells_AUC, file="cells_AUC.RData")
                       set.seed(123)
                       par(mfrow=c(1,1))
                       cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
                     })
        
      })
      
      output$cell_assign_score <- downloadHandler(
        filename = function()
          { paste("data-", Sys.Date(), ".csv", sep="")
          },
        content = function(file){

          write.csv(AUCell_table(),file)


        })




    output$cell_assign <- downloadHandler(
    filename = function(){'AUC.pdf'},
    content = function(file){
      pdf(file)
      print(

        AUCell()
      )
      dev.off()
    })



    output$Stouffer_sc<-renderTable({
      if(input$action1==0){
        return()
      }
      withProgress(message = 'Calculating Wilcoxon P-Value...', value = 0,
                   {         
      matric<-pb@assays$RNA@counts

      cells_sum <- Matrix::colSums(matric>=3)

      matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

      matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

      cells_sum<-Matrix::rowSums(t(matric))

      matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
      r_name<-rownames(matric)

      #Stouffer_CellEnrich()
      media_genes<-abc[[input$st_sc]]@geneIds
      media_genes
      same_genes<-intersect(r_name,media_genes)
      if(length(same_genes)==1)
        message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
      if(length(same_genes)==0)
        message="no common gene "
      if(length(same_genes)>1){
        matric= matric[same_genes,]####### didnot work why?

        #Further insuring have cells with atleast 10% expressed genes
        matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
        print(dim(matric))
        seurat_RNA_matric<-matric
        print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
        print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
        print(dim(seurat_RNA_matric))
        ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
        seurat_RNA_matric[seurat_RNA_matric==0]=1
        seurat_RNA_matric<-log2(seurat_RNA_matric)
        seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
        ###Stouffer score and then boxplot for each cell type
        stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
        print(sum(is.na(stouffer_score)))
        cell_types<-Idents(pb)[names(stouffer_score)]
        unique_cell_types=unique(cell_types)
        Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
        ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
          geom_boxplot(position=position_dodge(.2)) +
          geom_jitter(shape=16, position=position_jitter(.1)) +
          theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}

      p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
      colnames(p_val_stouffer_score)=unique_cell_types
      rownames(p_val_stouffer_score)=unique_cell_types
      for (i in 1:dim(p_val_stouffer_score)[1])
      {
        for (j in 1:dim(p_val_stouffer_score)[2])
        {
          print(paste(i,"_",j,sep=""))
          a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
          b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
          p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
        }
      }

      format(round(p_val_stouffer_score, 4))
})
    })


    stouffer_cellenrich<-reactive({
      matric<-pb@assays$RNA@counts
      
      cells_sum <- Matrix::colSums(matric>=3)
      
      matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]
      
      matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]
      
      cells_sum<-Matrix::rowSums(t(matric))
      
      matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
      r_name<-rownames(matric)
      
      #Stouffer_CellEnrich()
      media_genes<-abc[[input$st_sc]]@geneIds
      media_genes
      same_genes<-intersect(r_name,media_genes)
      if(length(same_genes)==1)
        message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
      if(length(same_genes)==0)
        message="no common gene "
      if(length(same_genes)>1){
        matric= matric[same_genes,]####### didnot work why?
        
        #Further insuring have cells with atleast 10% expressed genes
        matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
        print(dim(matric))
        seurat_RNA_matric<-matric
        print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
        print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
        print(dim(seurat_RNA_matric))
        ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
        seurat_RNA_matric[seurat_RNA_matric==0]=1
        seurat_RNA_matric<-log2(seurat_RNA_matric)
        seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
        ###Stouffer score and then boxplot for each cell type
        stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
        print(sum(is.na(stouffer_score)))
        cell_types<-Idents(pb)[names(stouffer_score)]
        unique_cell_types=unique(cell_types)
        Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
        ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
          geom_boxplot(position=position_dodge(.2)) +
          geom_jitter(shape=16, position=position_jitter(.1)) +
          theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
      
    })
    
  output$au_score<-renderPrint({
    if(input$au1==0){
      return()
    }
    withProgress(message = 'Cell Assignments...', value = 0,
                 {
    au_mat<-pb@assays$RNA@counts
    au_mat<-as.matrix(au_mat)
    par("mar")

    cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
    cells_rankings
    #save(cells_rankings, file="cells_rankings.RData")
    cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
    #save(cells_AUC, file="cells_AUC.RData")
    set.seed(123)
    par(mfrow=c(1,1))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC,  assign=TRUE)
    cells_assignment
})
  })
  # output$geneset_making<-renderPrint({
  #   if(input$cell_geneset==0){
  #     return()
  #   }
  #   feta_name=c("Chemosensory - Signature","Astrocyte68","Antigen presenting cell (RPS high)41","Adrenal gland inflammatory cell102","AT2 cell30","b_cell","b_cell_plasmocyte","Basal_cell","CB CD34+23","Chondrocyte99","Dendritic_cell ","Endothelial_cell  ","Endothelial cell (endothelial to mesenchymal transition)66","Endothelial cell (APC)8","Enterocyte_progenitor","Enterocyte","Epithelial cell (intermediated)60","Epithelial_cell","Erythroid_cell ","Erythroid progenitor cell (RP high)12","Fasciculata_cell ","fetal B cells","fetal Basophil_Mast","fetal CD1C+ DCs","fetal CEPs","fetal CLEC9A+ DCs","fetal Collecting duct lineage","fetal Connecting tubule lineage","fetal Distal tubule lineage","fetal EBMPs","fetal EEPs","fetal Enteroendocrine cells_intestine_pancreas","fetal Enteroendocrine cells_stomach","fetal Erythroblast","fetal ETDs","fetal HSCs","fetal HSPCs","fetal IL1B+ Microglia","fetal ILC 3","fetal Islet beta cells","fetal Islet delta cells","fetal LEC","fetal Loop of Henle lineage","fetal Macrophages","fetal Meg progenitors","fetal Meg","fetal Megakaryoblasts","fetal Microglia","fetal Nephron progenitor cell","fetal NK cells","fetal pDCs","fetal Perivascular macrophages","fetal Phagocytic macrophages","fetal Plasma cells","fetal Podocyte lineage","fetal Proximal tubule lineage","fetal PTPRC+ Microglia","fetal Pulmonary neuroendocrine cells","fetal Renal vesicle","fetal S100A9+ DCs","fetal T cells","fetal TMEM119+ Microglia","fetal TRAF1+ APCs","fetal Ureteric tip","fetal Ureteric trunk","fetalAntigen-presenting macrophages","fetalEndocardium","VEC-adrenal","VEC-brain","VEC-kidney","VEC-liver","VEC-lung","VEC-placenta","VEC-spleen","Fetal acinar cell71","Fetal chondrocyte43","Fetal endocrine cell88","Fetal enterocyte 15","Fetal epithelial progenitor1","Fetal fibroblast17","Fetal skeletal muscle cell64","Fetal_mesenchymal_progenitor","Fetal_neuron ","Fetal_stromal_cell","Fibroblast","Gastric chief cell50","Gastric endocrine cell84","Goblet_cell ","Hepatocyte_Endodermal cell16","hESC87","Immature sertoli cell (Pre-Sertoli cell)92","Intercalated cell56","Intermediated cell97","Kidney intercalated cell101","Loop of Henle25","M2 Macrophage47","Neutrophil (RPS high)40","Neutrophil","Macrophage","Mast cell86","Myeloid cell93","Oligodendrocyte28","Pancreas exocrine cell44","Primordial germ cell62","Proximal tubule progenitor83","Proliferating T cell52","Sinusoidal endothelial cell61","Stromal_cell","Smooth_muscle_cell","Stratified_epithelial_cell","T_cell","Thyroid follicular cell32","Ventricle cardiomyocyte74","Ureteric bud cell77")
  #   abc=list()
  #   for( i in feta_name)
  #   {
  #     fetal_name=read.csv(paste("./",i,".tsv", sep =""))
  #     fetal_name<-unique(fetal_name[,1])
  #     fetal_list=as.vector(fetal_name)
  # 
  #     abc[[i]]=GeneSet(fetal_list,setName=i)
  #   }
  #   gensets<-GeneSetCollection(abc)
  # })

  output$cell_rank<-renderPlot({
    if(input$action1==0){
      return()
    }
    withProgress(message = 'Calculating Cell Ranks...', value = 0,
                 {

                
    au_mat<-pb@assays$RNA@counts
    au_mat<-as.matrix(au_mat)
    par("mar")

  #  geneSet<-GeneSetCollection(list(gs1,gs2,gs3,gs4,gs5,gs6,gs7,gs8,gs9,gs61,gs10,gs11,gs62,gs12,gs13,gs14,gs15,gs16,gs17,gs18,gs19,gs20,gs60,gs21,gs22,gs23,gs24,gs25,gs26,gs27,gs28,gs29,gs30,gs31,gs32,gs33,gs34,gs35,gs36,gs37,gs38,gs39,gs40,gs41,gs42,gs43,gs44,gs45,gs46,gs47,gs48,gs49,gs50,gs51,gs52,gs53,gs54,gs55,gs56,gs57))


    cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
    cells_AUC <- AUCell_calcAUC(gensets, cells_rankings)
})

  })
  output$cell_auc<-renderPrint({
    if(input$cell_Auc==0){
      return()
    }
    withProgress(message = 'Calculating Cell AUC...', value = 0,
                 {        
    au_mat<-pb@assays$RNA@counts
    au_mat<-as.matrix(au_mat)
    #par("mar")

   # geneSet<-GeneSetCollection(list(gs1,gs2,gs3,gs4,gs5,gs6,gs7,gs8,gs9,gs61,gs10,gs11,gs62,gs12,gs13,gs14,gs15,gs16,gs17,gs18,gs19,gs20,gs60,gs21,gs22,gs23,gs24,gs25,gs26,gs27,gs28,gs29,gs30,gs31,gs32,gs33,gs34,gs35,gs36,gs37,gs38,gs39,gs40,gs41,gs42,gs43,gs44,gs45,gs46,gs47,gs48,gs49,gs50,gs51,gs52,gs53,gs54,gs55,gs56,gs57))


    cells_rankings <- AUCell_buildRankings(au_mat, nCores=1)

    cells_AUC <- AUCell_calcAUC(gensets, cells_rankings)
    cells_AUC
})

})
  output$hmap<-renderPlot({
    if(input$au4==0){
      return()
    }
    au_mat<-pb@assays$RNA@counts
    au_mat<-as.matrix(au_mat)
    par("mar")

   # geneSet<-GeneSetCollection(list(gs1,gs2,gs3,gs4,gs5,gs6,gs7,gs8,gs9,gs61,gs10,gs11,gs62,gs12,gs13,gs14,gs15,gs16,gs17,gs18,gs19,gs20,gs60,gs21,gs22,gs23,gs24,gs25,gs26,gs27,gs28,gs29,gs30,gs31,gs32,gs33,gs34,gs35,gs36,gs37,gs38,gs39,gs40,gs41,gs42,gs43,gs44,gs45,gs46,gs47,gs48,gs49,gs50,gs51,gs52,gs53,gs54,gs55,gs56,gs57))


    cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
    cells_rankings
  #  save(cells_rankings, file="cells_rankings.RData")
    cells_AUC <- AUCell_calcAUC(gensets, cells_rankings)
    #save(cells_AUC, file="cells_AUC.RData")
    set.seed(123)
    par(mfrow=c(1,1))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC,  assign=TRUE)
     cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
    assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
    colnames(assignmentTable)[2] <- "gensets"
    head(assignmentTable)
    assignmentMat <- table(assignmentTable[,"gensets"], assignmentTable[,"cell"])
    assignmentMat[,1:2]
    set.seed(123)
    miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]
    library(NMF)
    aheatmap(miniAssigMat, scale="none", color="pink")
  })


output$tsne<-renderPlot({
  if(input$action1==0){
    return()
  }
  pb <<- RunUMAP(pb, dims = 1:10)
  DimPlot(pb, reduction = "umap",pt.size = 3)

})

  output$c_tsne<-renderPlot({
    if(input$action1==0){
      return()
    }
    withProgress(message = 'Calculating Cell Signatures...', value = 0,{
                 
    au_mat<-pb@assays$RNA@counts
    au_mat<-as.matrix(au_mat)



    cells_rankings <- AUCell_buildRankings(au_mat, nCores=1, plotStats=TRUE)
    cells_rankings
    #save(cells_rankings, file="cells_rankings.RData")
    cells_AUC <- AUCell_calcAUC(abc[[input$aucell]], cells_rankings)
    cells_AUC
    set.seed(123)
    par(mfrow=c(1,1))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
    sumByGene <- apply(au_mat, 1, sum)
    exprMatSubset <- au_mat[which(sumByGene>0),]
    logMatrix <- log2(exprMatSubset+1)

    # umap_mat<-pb[["umap"]]@cell.embeddings
    umap_mat<<-as.data.frame(umap_mat)
    selectedThresholds <- getThresholdSelected(cells_assignment)

    for(geneSetName in names(selectedThresholds))
    {
      nBreaks <- 5 # Number of levels in the color palettes
      # Color palette for the cells that do not pass the threshold
      colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
      # Color palette for the cells that pass the threshold
      colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)

      # Split cells according to their AUC value for the gene set
      passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
      if(sum(passThreshold) >0 )
      {
        aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)

        # Assign cell color
        cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
                       setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))

        # Plot
        plot(umap_mat, main=geneSetName,
             sub="Pink/red cells pass the threshold",
             col=cellColor[rownames(umap_mat)], pch=16)
      }
    }
})
    
  })

  output$cell_umap <- downloadHandler(
    filename = function(){'umap.pdf'},
    content = function(file){
      pdf(file)
      print(
        # Plot

            plot(umap_mat, main=geneSetName,
                 sub="Pink/red cells pass the threshold",
                 col=cellColor[rownames(umap_mat)], pch=16)



      )
      dev.off()
    })
  output$cell_geneno <- downloadHandler(
    filename = function(){'combined.pdf'},
    content = function(file){
      pdf(file)
      print(
         plot(umap_mat, axes=FALSE, xlab="", ylab="",
       sub="Number of detected genes",
       col=cellColorNgenes[rownames(umap_mat)], pch=16)
  # Plot
      )

      dev.off()
    })


  output$hpa_stouffer<-renderPlot({
    if(input$action1==0){
      return()
    }
    withProgress(message = 'Generating Stouffer Plot...', value = 0,{
      
    matric<-pb@assays$RNA@counts

    cells_sum <- Matrix::colSums(matric>=3)

    matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

    matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

    cells_sum<-Matrix::rowSums(t(matric))

    matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
    r_name<-rownames(matric)
    media_genes <- pqr[[input$hpa11]]@geneIds
    same_genes<-intersect(r_name,media_genes)
    if(length(same_genes)==1)
      message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
    if(length(same_genes)==0)
      message="no common gene "
    if(length(same_genes)>1){
      matric= matric[same_genes,]####### didnot work why?

      #Further insuring have cells with atleast 10% expressed genes
      matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
      print(dim(matric))
      seurat_RNA_matric<-matric
      print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
      print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
      print(dim(seurat_RNA_matric))
      ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
      seurat_RNA_matric[seurat_RNA_matric==0]=1
      seurat_RNA_matric<-log2(seurat_RNA_matric)
      seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
      ###Stouffer score and then boxplot for each cell type
      stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
      print(sum(is.na(stouffer_score)))
      cell_types<-Idents(pb)[names(stouffer_score)]
      unique_cell_types=unique(cell_types)
      Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
      ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
        geom_boxplot(position=position_dodge(.2)) +
        geom_jitter(shape=16, position=position_jitter(.1)) +
        theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
  })
  })

  output$hpa_set2<-renderTable({
    if(input$action1==0){
      return()
    }
    withProgress(message = 'Calculating Wilcoxon P-Value...', value = 0,{
      
    matric<-pb@assays$RNA@counts

    cells_sum <- Matrix::colSums(matric>=3)

    matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

    matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

    cells_sum<-Matrix::rowSums(t(matric))

    matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
    r_name<-rownames(matric)
    media_genes <- pqr[[input$hpa11]]@geneIds
    same_genes<-intersect(r_name,media_genes)
    if(length(same_genes)==1)
      message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
    if(length(same_genes)==0)
      message="no common gene "
    if(length(same_genes)>1){
      matric= matric[same_genes,]####### didnot work why?

      #Further insuring have cells with atleast 10% expressed genes
      matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
      print(dim(matric))
      seurat_RNA_matric<-matric
      print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
      print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
      print(dim(seurat_RNA_matric))
      ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
      seurat_RNA_matric[seurat_RNA_matric==0]=1
      seurat_RNA_matric<-log2(seurat_RNA_matric)
      seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
      ###Stouffer score and then boxplot for each cell type
      stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
      print(sum(is.na(stouffer_score)))
      cell_types<-Idents(pb)[names(stouffer_score)]
      unique_cell_types=unique(cell_types)
      Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
      ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
        geom_boxplot(position=position_dodge(.2)) +
        geom_jitter(shape=16, position=position_jitter(.1)) +
        theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
    p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
    colnames(p_val_stouffer_score)=unique_cell_types
    rownames(p_val_stouffer_score)=unique_cell_types
    for (i in 1:dim(p_val_stouffer_score)[1])
    {
      for (j in 1:dim(p_val_stouffer_score)[2])
      {
        print(paste(i,"_",j,sep=""))
        a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
        b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
        p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
      }
    }
    format(round(p_val_stouffer_score, 4))
})
  })


 output$gtex_stouffer<-renderPlot({
   if(input$action1==0){
     return()
   }
   withProgress(message = 'Generating Stouffer Plot...', value = 0,{
     
   matric<-pb@assays$RNA@counts

   cells_sum <- Matrix::colSums(matric>=3)

   matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

   matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

   cells_sum<-Matrix::rowSums(t(matric))

   matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
   r_name<-rownames(matric)
   media_genes <- xyz[[input$gtex11]]@geneIds
   same_genes<-intersect(r_name,media_genes)
   if(length(same_genes)==1)
     message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
   if(length(same_genes)==0)
     message="no common gene "
   if(length(same_genes)>1){
     matric= matric[same_genes,]####### didnot work why?

     #Further insuring have cells with atleast 10% expressed genes
     matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
     print(dim(matric))
     seurat_RNA_matric<-matric
     print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
     print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
     print(dim(seurat_RNA_matric))
     ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
     seurat_RNA_matric[seurat_RNA_matric==0]=1
     seurat_RNA_matric<-log2(seurat_RNA_matric)
     seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
     ###Stouffer score and then boxplot for each cell type
     stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
     print(sum(is.na(stouffer_score)))
     cell_types<-Idents(pb)[names(stouffer_score)]
     unique_cell_types=unique(cell_types)
     Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
     ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
       geom_boxplot(position=position_dodge(.2)) +
       geom_jitter(shape=16, position=position_jitter(.1)) +
       theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
   })
})


 output$gtex_pstouffer<-renderTable({
   if(input$action1==0){
     return()
   }
   withProgress(message = 'Calculating Stouffer P-Value...', value = 0,{
     
   matric<-pb@assays$RNA@counts

   cells_sum <- Matrix::colSums(matric>=3)

   matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

   matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

   cells_sum<-Matrix::rowSums(t(matric))

   matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
   r_name<-rownames(matric)
   media_genes <- xyz[[input$gtex11]]@geneIds
   same_genes<-intersect(r_name,media_genes)
   if(length(same_genes)==1)
     message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
   if(length(same_genes)==0)
     message="no common gene "
   if(length(same_genes)>1){
     matric= matric[same_genes,]####### didnot work why?

     #Further insuring have cells with atleast 10% expressed genes
     matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
     print(dim(matric))
     seurat_RNA_matric<-matric
     print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
     print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
     print(dim(seurat_RNA_matric))
     ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
     seurat_RNA_matric[seurat_RNA_matric==0]=1
     seurat_RNA_matric<-log2(seurat_RNA_matric)
     seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
     ###Stouffer score and then boxplot for each cell type
     stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
     print(sum(is.na(stouffer_score)))
     cell_types<-Idents(pb)[names(stouffer_score)]
     unique_cell_types=unique(cell_types)
     Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
     ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
       geom_boxplot(position=position_dodge(.2)) +
       geom_jitter(shape=16, position=position_jitter(.1)) +
       theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}


   p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
   colnames(p_val_stouffer_score)=unique_cell_types
   rownames(p_val_stouffer_score)=unique_cell_types
   for (i in 1:dim(p_val_stouffer_score)[1])
   {
     for (j in 1:dim(p_val_stouffer_score)[2])
     {
       print(paste(i,"_",j,sep=""))
       a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
       b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
       p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
     }
   }
   format(round(p_val_stouffer_score, 4))
   })

 })

output$ex_cellrank<-renderPrint({
  if(input$action1==0){
    return()
  }


paste("The first step to calculate the enrichment of a signature is to create the rankings. These rankings are only an intermediate step to calculate the AUC. For each cell, the genes are ranked from highest to lowest value. The genes with the same expression value are shuffled. Therefore, genes with the expression '0' are randomly sorted at the end of the ranking.")
})
output$ex_AUC<-renderPrint({
  if(input$action1==0){
    return()
  }
paste("To determine whether the gene set is enriched at the top of the gene-ranking for each cell, AUCell uses the Area Under the Curve (AUC) of the recovery curve. In order to calculate the AUC, by default, only the top 5% of the genes in the ranking are used (i.e. checks whether the genes in the gene-set or signature are within the top 5%). This allows faster execution on bigger datasets and reduces the effect of the noise at the bottom of the ranking (e.g. where many genes might be tied at 0 counts). ")
    })
  output$ex_assign<-renderPrint({
    if(input$action1==0){
      return()
    }
paste("The AUC represents the proportion of expressed genes in the signature and their relative expression value compared to the other genes within the cell. We can use this property to explore the population of cells that are present in the dataset according to the expression of the gene-set. To ease the selection of an assignment threshold, this function adjusts the AUCs of each gene-set to several distributions and calculates possible thresholds: a) minimumDens (plot in Blue): Inflection point of the density curve. This is usually a good option for an ideal situation with bimodal distributions. To avoid false positives, by default this threshold will not be chosen if the second distribution is higher (i.e. the majority of cells have the gene-set active). b) L_k2 (plot in Red): Left distribution, after adjusting the AUC to a mixture of two distributions. The threshold is set to the right (prob: 1-(thrP/nCells)). c) R_k3 (plot in Pink): Right distribution, after adjusting the AUC to a mixture of three distributions. The threshold is set to the left (prob: thrP). d) Global_k1 (plot in Grey): global distribution (i.e. mean and standard deviations of all cells). The threshold is set to the right (prob: 1-(thrP/nCells))")
      })
output$ex_heat<-renderPrint({
  if(input$action1==0){
    return()
  }
paste("Extract these cells for all the gene-sets and transform it into a table of 0 and 1, where 1 represents the presence of that cell in a particular geneset and the heatmap is plotted based on that.")
  })
output$ex_umappp<-renderPrint({
  if(input$action1==0){
    return()
  }
paste("UMAP is colored based on the AUC scores. To highlight the cluster of cells that are more likely of the cell type according to the signatures, we split the cells into the cells that pass the assignment threshold (colored in shades of pink-red), and the cells that don't (colored in black-blue)")
  })
output$ex_uma<-renderPrint({
  if(input$action1==0){
    return()
  }
  paste("  The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. ")
})
output$ex_stouffer<-renderPrint({
  if(input$action1==0){
    return()
  }
paste("Stouffer's Z, a method closely related to Fisher's method, is based on Z-scores rather than p-values, allowing incorporation of study weights. It is named after the sociologist, Samuel A. Stouffer. ")
  })
output$n_gene<-renderPlot({
  if(input$ngene==0){
    return()
  }

   au_mat<-pb@assays$RNA@counts
  au_mat<-as.matrix(au_mat)

  pb$celltype<-Idents(pb)
  umap_mat<-pb[["umap"]]@cell.embeddings
  umap_mat<-as.data.frame(umap_mat)

  ###########################

  nGenesPerCell <- apply(au_mat, 2, function(x) sum(x>0))
  colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
  cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=FALSE,include.lowest=TRUE))], names(nGenesPerCell))
  plot(umap_mat, axes=FALSE, xlab="", ylab="",
       main="Number of Detected Genes",
       col=cellColorNgenes[rownames(umap_mat)], pch=16)




})
output$hpa_grn_1<- renderSimpleNetwork({
  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Network for HPA Datasets...', value = 0,{
    
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  total_tissue_gene<-gene_name_cluster
   user_genes<-gene_name_cluster


  # tissue_name=list("heart_human","placenta_human","placenta_mouse")

   hpa_cell=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")
   
   for( i in hpa_cell)
   {
     s<-pqr[[i]]@geneIds
     
     if(length(pqr[[i]]@geneIds)=="NULL")
       message=paste("null")
     
     # if(s@geneIds == "NULL")
     #   message=paste("no common genes")
     tissue_genes=s
     # if(length(tissue_genes)==0)
     #   message=paste("no common genes",sep="")
     genes_list=intersect(tissue_genes,user_genes)
     # total_tissue<-scan("tissue_genes",character())
     #total_tissue_gene<-scan("user_genes",character())
     x<-nrow(pb)
     if(length(genes_list)>1)
       
     {
       #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
       temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
       hpa_cell=rbind(hpa_cell,temp)
     }
   }
   


  #P_fold__df<-subset(P_fold__df,P_fold__df$logPvalue!= 0.000000e+00 && P_fold__df$adjusatedPvalue!=0.000000e+00)

  hpa_common<-as.data.frame(hpa_cell)
 # pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "HPA Dataset Tissues"

  final_hpa_grn$mor <- as.factor(final_hpa_grn$mor)

  data_1<-data_frame(
    from=final_hpa_grn$Target,
    to=final_hpa_grn$`HPA Dataset Tissues`
  )


 # p <- simpleNetwork(data_1, height="100px", width="100px")
  p<-simpleNetwork(data_1, height="200px", width="300px",linkDistance = 200, charge = -50, fontSize = 10, fontFamily = "serif",
                   linkColour = "#FFB6C1", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

  #par(mfrow=c(1,1))


  #plotm <- as.matrix(final_hpa_grn)
  # plotg <- graph_from_edgelist(rbind(plotm[,1:2],plotm[,2:3]), directed = T)
  # l <- layout_with_sugiyama(plotg, ceiling(match(V(plotg)$name, plotm)/nrow(plotm)),vgap = 20,hgap = 20)
  #nvv_col <- brewer.pal(9, "Oranges")[6]
  # hello<- plot(plotg, layout=-l$layout[,2:1],vertex.size = 30, vertex.color = "PINK",
  #    vertex.frame.color = "RED", vertex.label.color = "black", vertex.label.dist=1, vertex.label.cex = 1.5,main="Gene Regulatory Network")

p
})
  
})
output$hpa_grn<- renderSimpleNetwork({
  if(input$action1==0){
    return()
  }

  withProgress(message = 'Gene Regulatory Network for HPA Datasets...', value = 0,{
    

    marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
    
    gene_name_cluster<- marker_dataframe$gene
    # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
    gene_name_cluster
    total_tissue_gene<-gene_name_cluster
    user_genes<-gene_name_cluster
    
    
    # tissue_name=list("heart_human","placenta_human","placenta_mouse")
    
    hpa_cell=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")
    
    for( i in hpa_cell)
    {
      s<-pqr[[i]]@geneIds
      
      if(length(pqr[[i]]@geneIds)=="NULL")
        message=paste("null")
      
      # if(s@geneIds == "NULL")
      #   message=paste("no common genes")
      tissue_genes=s
      # if(length(tissue_genes)==0)
      #   message=paste("no common genes",sep="")
      genes_list=intersect(tissue_genes,user_genes)
      # total_tissue<-scan("tissue_genes",character())
      #total_tissue_gene<-scan("user_genes",character())
      x<-nrow(pb)
      if(length(genes_list)>1)
        
      {
        #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
        temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
        hpa_cell=rbind(hpa_cell,temp)
      }
    }
    
    


  #P_fold__df<-subset(P_fold__df,P_fold__df$logPvalue!= 0.000000e+00 && P_fold__df$adjusatedPvalue!=0.000000e+00)

  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "HPA Dataset Tissues"

  final_hpa_grn$mor <- as.factor(final_hpa_grn$mor)

data_1<-data_frame(
  from=final_hpa_grn$Target,
  to=final_hpa_grn$`HPA Type`,
)

data_2<-data_frame(
  from=final_hpa_grn$`Transcription factors`,
  to=final_hpa_grn$Target,
)
#p <- simpleNetwork(data_1, height="100px", width="100px")
q<-simpleNetwork(data_2, height="200px", width="300px",linkDistance = 200, charge = -30, fontSize = 10, fontFamily = "serif",
                 linkColour = "#90EE90", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

  #par(mfrow=c(1,1))

q
  #plotm <- as.matrix(final_hpa_grn)
 # plotg <- graph_from_edgelist(rbind(plotm[,1:2],plotm[,2:3]), directed = T)
 # l <- layout_with_sugiyama(plotg, ceiling(match(V(plotg)$name, plotm)/nrow(plotm)),vgap = 20,hgap = 20)
  #nvv_col <- brewer.pal(9, "Oranges")[6]
 # hello<- plot(plotg, layout=-l$layout[,2:1],vertex.size = 30, vertex.color = "PINK",
           #    vertex.frame.color = "RED", vertex.label.color = "black", vertex.label.dist=1, vertex.label.cex = 1.5,main="Gene Regulatory Network")

})

})
output$hpa_grn_table<-renderTable({
  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Table for HPA Datasets...', value = 0,{
    
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  total_tissue_gene<-gene_name_cluster

  user_genes<-gene_name_cluster

  hpa_cell=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")
  
  for( i in hpa_cell)
  {
    s<-pqr[[i]]@geneIds
    
    if(length(pqr[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  
  # tissue_name=list("heart_human","placenta_human","placenta_mouse")


  
  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "HPA Dataset Tissues"
  final_hpa_grn
  # tissue_name=list("heart_human","placenta_human","placenta_mouse")
})
})
output$gtex_grn_1<-renderSimpleNetwork({

  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Network for GTEx Datasets...', value = 0,{
    
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  user_genes<-gene_name_cluster
  
  hpa_cell = c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix-Uterine","Pituitary","Vagina","Pancreas")
  
  for( i in hpa_cell)
  {
    s<-xyz[[i]]@geneIds
    
    if(length(xyz[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }

  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_gtex_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V1")] <- "Target"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V2")] <- "GTEx Datasets Tissue"

  final_gtex_grn$mor <- as.factor(final_gtex_grn$mor)

  final_gtex_grn


  data_1<-data_frame(
    from=final_gtex_grn$`Transcription factors`,
    to=final_gtex_grn$Target
  )


  p<-simpleNetwork(data_1, height="200px", width="300px",linkDistance = 200, charge = -50, fontSize = 10, fontFamily = "serif",
                   linkColour = "#90EE90", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)
p
})
})

output$gtex_grn<-renderSimpleNetwork({

  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Network for GTEx Datasets...', value = 0,{
    
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  user_genes<-gene_name_cluster
  
  hpa_cell = c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix-Uterine","Pituitary","Vagina","Pancreas")
  
  for( i in hpa_cell)
  {
    s<-xyz[[i]]@geneIds
    
    if(length(xyz[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  
  # tissue_name=list("heart_human","placenta_human","placenta_mouse")


  hpa_common<-as.data.frame(hpa_cell)
 # pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_gtex_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V1")] <- "Target"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V2")] <- "GTEx Dataset Tissues"

  final_gtex_grn$mor <- as.factor(final_gtex_grn$mor)

  final_gtex_grn


  data_1<-data_frame(
    from=final_gtex_grn$Target,
    to=final_gtex_grn$`GTEx Dataset Tissues`
  )

  
  ######
  

  # p <- simpleNetwork(data_1, height="100px", width="100px")
  p<-simpleNetwork(data_1, height="200px", width="300px",linkDistance = 200, charge = -50, fontSize = 10, fontFamily = "serif",
                   linkColour = "	#FFB6C1", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

p
})
})
output$gtex_grn_table<-renderTable({

  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Table for GTEx Datasets...', value = 0,{
    
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  user_genes<-gene_name_cluster
  
  hpa_cell = c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix-Uterine","Pituitary","Vagina","Pancreas")
  
  for( i in hpa_cell)
  {
    s<-xyz[[i]]@geneIds
    
    if(length(xyz[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  

  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_gtex_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V1")] <- "Target"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V2")] <- "GTEx Dataset Tissues "

  final_gtex_grn$mor <- as.factor(final_gtex_grn$mor)

  final_gtex_grn
})
  
})

output$cell_grn_1<-renderSimpleNetwork({
  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Network...', value = 0,{
    marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
    
    gene_name_cluster<- marker_dataframe$gene
    # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
    gene_name_cluster
    total_tissue_gene<-gene_name_cluster
  ###########
    user_genes<-gene_name_cluster
    hpa_cell =c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")
    for( i in hpa_cell)
    {
      s<-abc[[i]]@geneIds
      
      if(length(abc[[i]]@geneIds)=="NULL")
        message=paste("null")
      
      # if(s@geneIds == "NULL")
      #   message=paste("no common genes")
      tissue_genes=s
      # if(length(tissue_genes)==0)
      #   message=paste("no common genes",sep="")
      genes_list=intersect(tissue_genes,user_genes)
      # total_tissue<-scan("tissue_genes",character())
      #total_tissue_gene<-scan("user_genes",character())
      x<-nrow(pb)
      if(length(genes_list)>1)
        
      {
        #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
        temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
        hpa_cell=rbind(hpa_cell,temp)
      }
    }
    
  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_cell_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V1")] <- "Target"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V2")] <- "Cell Enrich Signature"

  final_cell_grn$mor <- as.factor(final_cell_grn$mor)


  data_2<-data_frame(
    from=final_cell_grn$`Transcription factors`,
    to=final_cell_grn$Target
  )


  #p <- simpleNetwork(data_1, height="100px", width="100px")
  q<-simpleNetwork(data_2, height="200px", width="300px",linkDistance = 200, charge = -30, fontSize = 10, fontFamily = "serif",
                   linkColour = "#90EE90", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

  ######


q
})

})






output$cell_grn<-renderSimpleNetwork({
  if(input$action1==0){
    return()
  }
  withProgress(message = 'Gene Regulatory Network...', value = 0,{
    
    marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
    
    gene_name_cluster<- marker_dataframe$gene
    # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
    gene_name_cluster
    total_tissue_gene<-gene_name_cluster
    user_genes<-gene_name_cluster
  ###########
    hpa_cell =c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")
    for( i in hpa_cell)
    {
      s<-abc[[i]]@geneIds
      
      if(length(abc[[i]]@geneIds)=="NULL")
        message=paste("null")
      
      # if(s@geneIds == "NULL")
      #   message=paste("no common genes")
      tissue_genes=s
      # if(length(tissue_genes)==0)
      #   message=paste("no common genes",sep="")
      genes_list=intersect(tissue_genes,user_genes)
      # total_tissue<-scan("tissue_genes",character())
      #total_tissue_gene<-scan("user_genes",character())
      x<-nrow(pb)
      if(length(genes_list)>1)
        
      {
        #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
        temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
        hpa_cell=rbind(hpa_cell,temp)
      }
    }
    

  hpa_common<-as.data.frame(hpa_cell)
 # pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_cell_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V1")] <- "Target"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V2")] <- "Cell Enrich Signature"

  final_cell_grn$mor <- as.factor(final_cell_grn$mor)


  data_1<-data_frame(
    from=final_cell_grn$Target,
    to=final_cell_grn$`Cell Enrich Signature`
  )


  #p <- simpleNetwork(data_1, height="100px", width="100px")
  q<-simpleNetwork(data_1, height="200px", width="300px",linkDistance = 200, charge = -30, fontSize = 10, fontFamily = "serif",
                   linkColour = "#FFB6C1", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

  ######

q

})

})



















cell_grn_d<-reactive({
  tissue_name=c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")



  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster

  total_tissue_gene<-gene_name_cluster
  user_genes<-gene_name_cluster


  # tissue_name=list("heart_human","placenta_human","placenta_mouse")

  hpa_cell=c()

  for( i in tissue_name)
  {
    tissue_genes=read.csv(paste("./",i,".tsv", sep =""))
    genes_list=intersect(tissue_genes[,1],user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)

    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)

      #f1<-(length_overlap/length_tissue_genes)
      #f2<- (tissue_genes_length/20400)
      #fold_change<-f1/f2 incorrect number of subscripts on matrix
    }
  }

  # tissue_name=list("heart_human","placenta_human","placenta_mouse")


  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_cell_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V1")] <- "Target"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V2")] <- "Cell Enrich Signature"

  final_cell_grn$mor <- as.factor(final_cell_grn$mor)
  ######################

})


output$cell_grn_table<-renderTable({
  if(input$action1==0){
    return()
  }
  withProgress(message = 'Generating Gene Regulatory Table... ', value = 0,{
    
    marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
    
    gene_name_cluster<- marker_dataframe$gene
    # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
    gene_name_cluster
    total_tissue_gene<-gene_name_cluster
    user_genes<-gene_name_cluster          
    hpa_cell =c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")
    for( i in hpa_cell)
    {
      s<-abc[[i]]@geneIds
      
      if(length(abc[[i]]@geneIds)=="NULL")
        message=paste("null")
      
      # if(s@geneIds == "NULL")
      #   message=paste("no common genes")
      tissue_genes=s
      # if(length(tissue_genes)==0)
      #   message=paste("no common genes",sep="")
      genes_list=intersect(tissue_genes,user_genes)
      # total_tissue<-scan("tissue_genes",character())
      #total_tissue_gene<-scan("user_genes",character())
      x<-nrow(pb)
      if(length(genes_list)>1)
        
      {
        #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
        temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
        hpa_cell=rbind(hpa_cell,temp)
      }
    }
    
  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "Cell Enrich Signature"
  final_hpa_grn
  })# tissue_name=list("heart_human","placenta_human","placenta_mouse")

})



output$cell_stoudownload <- downloadHandler(
  filename = function(){'Stouffer CellEnrich.pdf'},
  content = function(file){
    pdf(file)
    print(


stouffer_cellenrich()

    )
    dev.off()
  })


cell_grn_d<-reactive({
  tissue_name=c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")


  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster


  total_tissue_gene<-gene_name_cluster
  user_genes<-gene_name_cluster


  # tissue_name=list("heart_human","placenta_human","placenta_mouse")

  hpa_cell=c()

  for( i in tissue_name)
  {
    tissue_genes=read.csv(paste("R/",i,".tsv", sep =""))
    genes_list=intersect(tissue_genes[,1],user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)

    {
     # write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)

      #f1<-(length_overlap/length_tissue_genes)
      #f2<- (tissue_genes_length/20400)
      #fold_change<-f1/f2 incorrect number of subscripts on matrix
    }
  }

  # tissue_name=list("heart_human","placenta_human","placenta_mouse")


  hpa_common<-as.data.frame(hpa_cell)
 # pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_cell_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V1")] <- "Target"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_cell_grn)[which(names(final_cell_grn) == "V2")] <- "Cell Enrich Signature"

  final_cell_grn$mor <- as.factor(final_cell_grn$mor)






  plotm <- as.matrix(final_cell_grn)
  plotg <- graph_from_edgelist(rbind(plotm[,1:2],plotm[,2:3]), directed = T)
  l <- layout_with_sugiyama(plotg, ceiling(match(V(plotg)$name, plotm)/nrow(plotm)),vgap = 20,hgap = 20)
  #nvv_col <- brewer.pal(9, "Oranges")[6]
  hello<- plot(plotg, layout=-l$layout[,2:1],vertex.size = 30, vertex.color = "PINK",
               vertex.frame.color = "RED", vertex.label.color = "black", vertex.label.dist=1, vertex.label.cex = 1.5,main="Gene Regulatory Network")
  ##############################


})



output$cell_grndownload <- downloadHandler(
  filename = function(){'GRN.pdf'},
  content = function(file){
    pdf(file)
    print(

      cell_grn_d()


    )
    dev.off()
  })



hpa_grn_d<-reactive({

  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  total_tissue_gene<-gene_name_cluster

  hpa_cell=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")
  
  for( i in hpa_cell)
  {
    s<-pqr[[i]]@geneIds
    
    if(length(pqr[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  
#P_fold__df<-subset(P_fold__df,P_fold__df$logPvalue!= 0.000000e+00 && P_fold__df$adjusatedPvalue!=0.000000e+00)

hpa_common<-as.data.frame(hpa_cell)
#pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

##########************************* read dorothea_regulon file
dorothea_regulon_human <-doro

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores
pb <- run_viper(pb, regulon)

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = pb) <- "dorothea"
pb <- ScaleData(pb)
pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                             logfc.threshold = input$log_fc , verbose = FALSE)

DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

## We transform Viper scores, scaled by seurat, into a data frame to better

## handling the results
viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(pb)),
                            cell_type = as.character(Idents(pb)),
                            stringsAsFactors = FALSE)


################## my addition ###################################
viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>%
  group_by(variable, cell_type) %>%
  summarise(avg = mean(value),
            std = sd(value))
## We select the 20 most variable TFs. (20*9 populations = 180)

highly_variable_tfs <- summarized_viper_scores %>%
  group_by(variable) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(variable)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "variable") %>%
  dplyr::select(-std) %>%
  spread(variable, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

transpose_summarise_df<-t(summarized_viper_scores_df)

merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "HPA Dataset Tissues"

final_hpa_grn$mor <- as.factor(final_hpa_grn$mor)

data_1<-data_frame(
  from=final_hpa_grn$Target,
  to=final_hpa_grn$`HPA Type`,
)


# p <- simpleNetwork(data_1, height="100px", width="100px")
p<-simpleNetwork(data_1, height="200px", width="300px",linkDistance = 200, charge = -50, fontSize = 10, fontFamily = "serif",
                 linkColour = "	#FF69B4", nodeColour = "#FFFFFF", opacity = 1, zoom = T	)

#par(mfrow=c(1,1))






})

gtex_d<-reactive({

  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
  
  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  user_genes<-gene_name_cluster
  
  hpa_cell = c("Adrenal Gland","Adipose Tissue","Brain","Colon","Esophagus","Fallopian Tube","Bladder","Heart", "Muscle","Kidney","Liver","Lung","Ovary","Prostate","Salivary Gland","Skin","Small Intestine","Spleen","Stomach","Thyroid","Breast","Nerve","Uterus","Cervix-Uterine","Pituitary","Vagina","Pancreas")
  
  for( i in hpa_cell)
  {
    s<-xyz[[i]]@geneIds
    
    if(length(xyz[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  
  
  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold
  
  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro
  
  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))
  
  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)
  
  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)
  
  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")
  
  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)
  
  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()
  
  ## We transform Viper scores, scaled by seurat, into a data frame to better
  
  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()
  
  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)
  
  
  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)
  
  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)
  
  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)
  
  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  transpose_summarise_df<-t(summarized_viper_scores_df)
  
  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_gtex_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V1")] <- "Target"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_gtex_grn)[which(names(final_gtex_grn) == "V2")] <- "GTEx Tissue "
  
  final_gtex_grn$mor <- as.factor(final_gtex_grn$mor)
  
  final_gtex_grn
})










cell_grn_dt<-reactive({
  marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )
  
  gene_name_cluster<- marker_dataframe$gene
  # write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
  gene_name_cluster
  total_tissue_gene<-gene_name_cluster
  user_genes<-gene_name_cluster
  hpa_cell =c("Ectodermal cell differentiation","Ectodermal development","Embryonic stem cell","Adult B cell","Adult B cell plasmocyte","Fetal fibroblast","Endodermal cell differentiation","Endodermal development","Fibroblast migration","Fibroblast proliferation","Mesoderm development","Mesodermal cell differentiation","Adult astrocyte","Adult antigen presenting cell (RPS high)","Adult adrenal gland inflammatory cell","Adult monocyte","Adult AT2 cell","Adult basal cell","Adult CB CD34+23","Developing heart atrial cardiomyocyte","Developing heart ventricular cardiomyocyte","Adult chondrocyte","Adult chemosensory","Adult dendritic cell","Adult endothelial cell","Adult endothelial cell (endothelial to mesenchymal transition)","Adult endothelial cell (APC)","Adult enterocyte progenitor","Adult enterocyte","Adult epithelial cell (intermediate)","Adult epithelial cell","Adult erythroid cell","Adult erythroid progenitor cell (RP high)","Adult fasciculata cell","Fetal B cells","Fetal basophil mast","Fetal CD1C+ DCs","Fetal CEPs","Fetal CLEC9A+ DCs","Fetal collecting duct lineage","Fetal connecting tubule lineage","Fetal distal tubule lineage","Fetal EBMPs","Fetal EEPs","Fetal enteroendocrine cells intestine pancreas","Fetal enteroendocrine cells stomach","Fetal erythroblast","Fetal ETDs","Fetal HSCs","Fetal HSPCs","Fetal IL1B+ microglia","Fetal ILC3","Fetal islet beta cells","Fetal islet delta cells","Fetal LEC","Fetal loop of Henle lineage","Fetal macrophages","Fetal meg progenitors","Fetal meg","Fetal megakaryoblasts","Fetal microglia","Fetal nephron progenitor cell","Fetal NK cells","Fetal pDCs","Fetal perivascular macrophages","Fetal phagocytic macrophages","Fetal plasma cells","Fetal podocyte lineage","Fetal proximal tubule lineage","Fetal PTPRC+ microglia","Fetal pulmonary neuroendocrine cells","Fetal renal vesicle","Fetal S100A9+ DCs","Fetal T cells","Fetal TMEM119+ microglia","Fetal TRAF1+ APCs","Fetal ureteric tip","Fetal ureteric trunk","Fetal antigen presenting macrophages","Fetal endocardium","Fetal adrenal","Fetal brain","Fetal kidney","Fetal liver","Fetal lung","Fetal placenta","Fetal spleen","Fetal acinar cell","Fetal chondrocyte","Fetal endocrine cell","Fetal enterocyte","Fetal epithelial progenitor","Fetal skeletal muscle cell","Fetal mesenchymal progenitor","Fetal neuron","Fetal stromal cell","Adult fibroblast","Adult gastric chief cell","Adult gastric endocrine cell","Adult goblet cell","Adult hepatocyte endodermal cell","Adult hESC87","Adult immature sertoli cell (pre-sertoli cell)","Adult intercalated cell","Adult intermediated cell","Adult kidney intercalated cell","Adult loop of Henle","Adult M2 macrophage","Adult neutrophil (RPS high)","Adult neutrophil","Adult macrophage","Adult mast cell","Adult myeloid cell","Adult oligodendrocyte","Adult pancreas exocrine cell","Adult primordial germ cell","Adult proximal tubule progenitor","Adult proliferating T cell","Adult sinusoidal endothelial cell","Adult stromal cell","Adult smooth muscle cell","Adult stratified epithelial cell","Adult T cell","Adult thyroid follicular cell","Adult ventricle cardiomyocyte","Adult ureteric bud cell")
  for( i in hpa_cell)
  {
    s<-abc[[i]]@geneIds
    
    if(length(abc[[i]]@geneIds)=="NULL")
      message=paste("null")
    
    # if(s@geneIds == "NULL")
    #   message=paste("no common genes")
    tissue_genes=s
    # if(length(tissue_genes)==0)
    #   message=paste("no common genes",sep="")
    genes_list=intersect(tissue_genes,user_genes)
    # total_tissue<-scan("tissue_genes",character())
    #total_tissue_gene<-scan("user_genes",character())
    x<-nrow(pb)
    if(length(genes_list)>1)
      
    {
      #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
      temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
      hpa_cell=rbind(hpa_cell,temp)
    }
  }
  
  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc , verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "Cell Enrich Signature"
  final_hpa_grn
  # tissue_name=list("


})
hpa_grn_dt<-reactive({ marker_dataframe <- subset(pb.markers,pb.markers$cluster == input$grn_obs )

gene_name_cluster<- marker_dataframe$gene
# write.table(gene_name_cluster,"gene_name_cluster.txt",row.names=F)
gene_name_cluster
  total_tissue_gene<-gene_name_cluster
    user_genes<-gene_name_cluster


  # tissue_name=list("heart_human","placenta_human","placenta_mouse")


    hpa_cell=c("Lymph Node-HPA",	"Tonsil-HPA",	"Appendix-HPA",	"Spleen-HPA",	"Bone Marrow-HPA"	,"Esophagus-HPA"	,"Skin-HPA",	"Colon-HPA",	"Rectum-HPA",	"Duodenum-HPA",	"Small Intestine-HPA",	"Stomach-HPA",	"Adipose Tissue-HPA",	"Lung-HPA",	"Placenta-HPA",	"Gallbladder-HPA",	"Urinary Bladder-HPA","Endometrium-HPA",	"Smooth Muscle-HPA","Fallopian Tube-HPA",	"Thyroid Gland-HPA",	"Ovary-HPA","Prostate",	"Kidney-HPA",	"Adrenal Gland-HPA",	"Brain-HPA",	"Salivary Gland-HPA","Pancreas-HPA",	"Skeletal Muscle-HPA",	"Liver-HPA",	"Heart Muscle-HPA")
    
    for( i in hpa_cell)
    {
      s<-pqr[[i]]@geneIds
      
      if(length(pqr[[i]]@geneIds)=="NULL")
        message=paste("null")
      
      # if(s@geneIds == "NULL")
      #   message=paste("no common genes")
      tissue_genes=s
      # if(length(tissue_genes)==0)
      #   message=paste("no common genes",sep="")
      genes_list=intersect(tissue_genes,user_genes)
      # total_tissue<-scan("tissue_genes",character())
      #total_tissue_gene<-scan("user_genes",character())
      x<-nrow(pb)
      if(length(genes_list)>1)
        
      {
        #write.csv(genes_list,paste(i,"_common_genes.csv",sep=""))
        temp=matrix(c(genes_list,rep(i,length(genes_list))),ncol=2)
        hpa_cell=rbind(hpa_cell,temp)
      }
    }
    


  #P_fold__df<-subset(P_fold__df,P_fold__df$logPvalue!= 0.000000e+00 && P_fold__df$adjusatedPvalue!=0.000000e+00)

  hpa_common<-as.data.frame(hpa_cell)
  #pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = input$log_fc)### have to make chnage in threshold

  ##########************************* read dorothea_regulon file
  dorothea_regulon_human <-doro

  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))

  ## We compute Viper Scores
  pb <- run_viper(pb, regulon)

  ## We compute the Nearest Neighbours to perform cluster
  DefaultAssay(object = pb) <- "dorothea"
  pb <- ScaleData(pb)
  pb <- RunPCA(pb, features = rownames(pb), verbose = FALSE)
  pb <- FindNeighbors(pb, dims = 1:10, verbose = FALSE)
  pb <- FindClusters(pb, resolution = 0.5, verbose = FALSE)

  pb <- RunUMAP(pb, dims = 1:10, umap.method = "uwot", metric = "cosine")

  pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = input$log_fc, verbose = FALSE)

  DimPlot(pb, reduction = "umap", label = FALSE, pt.size = 0.5) + NoLegend()

  ## We transform Viper scores, scaled by seurat, into a data frame to better

  ## handling the results
  viper_scores_df <- GetAssayData(pb, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()

  ## We create a data frame containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(pb)),
                              cell_type = as.character(Idents(pb)),
                              stringsAsFactors = FALSE)


  ################## my addition ###################################
  viper_scores_clusters <- viper_scores_df  %>% data.frame() %>% rownames_to_column("cell")
  viper_scores_clusters<- melt(viper_scores_clusters, id= "cell")
  viper_scores_clusters<- viper_scores_clusters %>% inner_join(CellsClusters)

  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(variable, cell_type) %>%
    summarise(avg = mean(value),
              std = sd(value))
  ## We select the 20 most variable TFs. (20*9 populations = 180)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(variable) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(180, var) %>%
    distinct(variable)

  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "variable") %>%
    dplyr::select(-std) %>%
    spread(variable, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  transpose_summarise_df<-t(summarized_viper_scores_df)

  merged.data <- merge(transpose_summarise_df, regulon,by.x="row.names", by.y="tf")
  final_data_plot<- merged.data  %>% dplyr::select(Row.names,confidence,target,mor)
  final_hpa_grn<-merge(hpa_common,final_data_plot,by.x="V1",by.y="target")
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V1")] <- "Target"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "Row.names")] <- "Transcription factors"
  colnames(final_hpa_grn)[which(names(final_hpa_grn) == "V2")] <- "HPA Dataset Tissues"

  final_hpa_grn$mor <- as.factor(final_hpa_grn$mor)

  final_hpa_grn


})





output$hpa_grndownload <- downloadHandler(
  filename = function(){'GRN.pdf'},
  content = function(file){
    pdf(file)
    print(
      hpa_grn_d()



    )
    dev.off()
  })
output$gtex_grndownload <- downloadHandler(
  filename = function(){'GRN.pdf'},
  content = function(file){
    pdf(file)
    print(

      gtex_grn_d()


    )
    dev.off()
  })

output$report <- downloadHandler(
  # For PDF output, change this to "report.pdf"
  filename = "report.html",
  content = function(file) {
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)

    # Set up parameters to pass to Rmd document
    params <- list(dataTables)

    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv())
    )
  }
)

stouffer_hpa<-reactive({
  matric<-pb@assays$RNA@counts

  cells_sum <- Matrix::colSums(matric>=3)

  matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

  matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

  cells_sum<-Matrix::rowSums(t(matric))

  matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
  r_name<-rownames(matric)
  media_genes <- pqr[[input$hpa11]]@geneIds
  same_genes<-intersect(r_name,media_genes)
  if(length(same_genes)==1)
    message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
  if(length(same_genes)==0)
    message="no common gene "
  if(length(same_genes)>1){
    matric= matric[same_genes,]####### didnot work why?

    #Further insuring have cells with atleast 10% expressed genes
    matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
    print(dim(matric))
    seurat_RNA_matric<-matric
    print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
    print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
    print(dim(seurat_RNA_matric))
    ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
    seurat_RNA_matric[seurat_RNA_matric==0]=1
    seurat_RNA_matric<-log2(seurat_RNA_matric)
    seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
    ###Stouffer score and then boxplot for each cell type
    stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
    print(sum(is.na(stouffer_score)))
    cell_types<-Idents(pb)[names(stouffer_score)]
    unique_cell_types=unique(cell_types)
    Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
    ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
      geom_boxplot(position=position_dodge(.2)) +
      geom_jitter(shape=16, position=position_jitter(.1)) +
      theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}



})

stouffer_gtex<-reactive({
  matric<-pb@assays$RNA@counts

  cells_sum <- Matrix::colSums(matric>=3)

  matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

  matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

  cells_sum<-Matrix::rowSums(t(matric))

  matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
  r_name<-rownames(matric)
  media_genes <- xyz[[input$gtex11]]@geneIds
  same_genes<-intersect(r_name,media_genes)
  if(length(same_genes)==1)
    message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
  if(length(same_genes)==0)
    message="no common gene "
  if(length(same_genes)>1){
    matric= matric[same_genes,]####### didnot work why?

    #Further insuring have cells with atleast 10% expressed genes
    matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
    print(dim(matric))
    seurat_RNA_matric<-matric
    print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
    print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
    print(dim(seurat_RNA_matric))
    ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
    seurat_RNA_matric[seurat_RNA_matric==0]=1
    seurat_RNA_matric<-log2(seurat_RNA_matric)
    seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
    ###Stouffer score and then boxplot for each cell type
    stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
    print(sum(is.na(stouffer_score)))
    cell_types<-Idents(pb)[names(stouffer_score)]
    unique_cell_types=unique(cell_types)
    Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
    ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
      geom_boxplot(position=position_dodge(.2)) +
      geom_jitter(shape=16, position=position_jitter(.1)) +
      theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
})




stouffer_tab_gtex<-reactive({

  matric<-pb@assays$RNA@counts

  cells_sum <- Matrix::colSums(matric>=3)

  matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

  matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

  cells_sum<-Matrix::rowSums(t(matric))

  matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
  r_name<-rownames(matric)
  media_genes <- xyz[[input$gtex11]]@geneIds
  same_genes<-intersect(r_name,media_genes)
  if(length(same_genes)==1)
    message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
  if(length(same_genes)==0)
    message="no common gene "
  if(length(same_genes)>1){
    matric= matric[same_genes,]####### didnot work why?

    #Further insuring have cells with atleast 10% expressed genes
    matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
    print(dim(matric))
    seurat_RNA_matric<-matric
    print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
    print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
    print(dim(seurat_RNA_matric))
    ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
    seurat_RNA_matric[seurat_RNA_matric==0]=1
    seurat_RNA_matric<-log2(seurat_RNA_matric)
    seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
    ###Stouffer score and then boxplot for each cell type
    stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
    print(sum(is.na(stouffer_score)))
    cell_types<-Idents(pb)[names(stouffer_score)]
    unique_cell_types=unique(cell_types)
    Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
    ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
      geom_boxplot(position=position_dodge(.2)) +
      geom_jitter(shape=16, position=position_jitter(.1)) +
      theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}


  p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
  colnames(p_val_stouffer_score)=unique_cell_types
  rownames(p_val_stouffer_score)=unique_cell_types
  for (i in 1:dim(p_val_stouffer_score)[1])
  {
    for (j in 1:dim(p_val_stouffer_score)[2])
    {
      print(paste(i,"_",j,sep=""))
      a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
      b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
      p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
    }
  }
  format(round(p_val_stouffer_score, 4))

})



output$download_hpa_gt <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv( hpa_grn_dt(), file)

  })



output$stouffer_table_gtex <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(stouffer_tab_gtex(), file)

  })

output$stouffer_table_hpa <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(stouffer_tab_hpa(), file)

  })
output$stouffer_table_cell <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(stouffer_tab_cell(), file)

  })
stouffer_tab_cell<-reactive({
  matric<-pb@assays$RNA@counts

  cells_sum <- Matrix::colSums(matric>=3)

  matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

  matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

  cells_sum<-Matrix::rowSums(t(matric))

  matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
  r_name<-rownames(matric)

  #Stouffer_CellEnrich()
  media_genes<-abc[[input$st_sc]]@geneIds
  media_genes
  same_genes<-intersect(r_name,media_genes)
  if(length(same_genes)==1)
    message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
  if(length(same_genes)==0)
    message="no common gene "
  if(length(same_genes)>1){
    matric= matric[same_genes,]####### didnot work why?

    #Further insuring have cells with atleast 10% expressed genes
    matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
    print(dim(matric))
    seurat_RNA_matric<-matric
    print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
    print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
    print(dim(seurat_RNA_matric))
    ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
    seurat_RNA_matric[seurat_RNA_matric==0]=1
    seurat_RNA_matric<-log2(seurat_RNA_matric)
    seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
    ###Stouffer score and then boxplot for each cell type
    stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
    print(sum(is.na(stouffer_score)))
    cell_types<-Idents(pb)[names(stouffer_score)]
    unique_cell_types=unique(cell_types)
    Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
    ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
      geom_boxplot(position=position_dodge(.2)) +
      geom_jitter(shape=16, position=position_jitter(.1)) +
      theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}

  p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
  colnames(p_val_stouffer_score)=unique_cell_types
  rownames(p_val_stouffer_score)=unique_cell_types
  for (i in 1:dim(p_val_stouffer_score)[1])
  {
    for (j in 1:dim(p_val_stouffer_score)[2])
    {
      print(paste(i,"_",j,sep=""))
      a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
      b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
      p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
    }
  }

  format(round(p_val_stouffer_score, 4))
})


stouffer_tab_hpa<-reactive({
  matric<-pb@assays$RNA@counts

  cells_sum <- Matrix::colSums(matric>=3)

  matric<-matric[,intersect(which(cells_sum>=stats::quantile(cells_sum,probs = 0.001)), which(cells_sum<=stats::quantile(cells_sum,probs = 1)))]

  matric<-matric[which(Matrix::rowSums(matric > 2) > 3),]

  cells_sum<-Matrix::rowSums(t(matric))

  matric<-Matrix::t(t(matric)/(cells_sum/stats::median(cells_sum)))
  r_name<-rownames(matric)
  media_genes <- pqr[[input$hpa11]]@geneIds
  same_genes<-intersect(r_name,media_genes)
  if(length(same_genes)==1)
    message=paste("There is only one common gene ",same_genes," so cannot compute stouffer score ",sep="")
  if(length(same_genes)==0)
    message="no common gene "
  if(length(same_genes)>1){
    matric= matric[same_genes,]####### didnot work why?

    #Further insuring have cells with atleast 10% expressed genes
    matric<-matric[,Matrix::colSums(matric>0)>(dim(matric)[1]/10)]
    print(dim(matric))
    seurat_RNA_matric<-matric
    print(sum(apply(seurat_RNA_matric,1,function(x) sum(x)==0)))
    print(sum(apply(seurat_RNA_matric,2,function(x) sum(x)==0)))
    print(dim(seurat_RNA_matric))
    ### First adding 1 (adding 1 to only zero, not only for zeros) then log2 then zscore
    seurat_RNA_matric[seurat_RNA_matric==0]=1
    seurat_RNA_matric<-log2(seurat_RNA_matric)
    seurat_RNA_matric<-t(apply(seurat_RNA_matric,1,function(x) x*mean(x)/sd(x)))
    ###Stouffer score and then boxplot for each cell type
    stouffer_score<- apply(seurat_RNA_matric,2,function(x) sum(x)/sqrt(length(x)))
    print(sum(is.na(stouffer_score)))
    cell_types<-Idents(pb)[names(stouffer_score)]
    unique_cell_types=unique(cell_types)
    Stouffer_score_df<-data.frame("Stouffer_score"=stouffer_score, "cell_types"=cell_types)
    ggplot(Stouffer_score_df, aes(x=cell_types, y=Stouffer_score, fill=cell_types)) +
      geom_boxplot(position=position_dodge(.2)) +
      geom_jitter(shape=16, position=position_jitter(.1)) +
      theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))}
  p_val_stouffer_score=matrix(0,nrow=length(unique_cell_types),ncol=length(unique_cell_types))
  colnames(p_val_stouffer_score)=unique_cell_types
  rownames(p_val_stouffer_score)=unique_cell_types
  for (i in 1:dim(p_val_stouffer_score)[1])
  {
    for (j in 1:dim(p_val_stouffer_score)[2])
    {
      print(paste(i,"_",j,sep=""))
      a=stouffer_score[names(cell_types)[cell_types==unique_cell_types[i]]]
      b=stouffer_score[names(cell_types)[cell_types==unique_cell_types[j]]]
      p_val_stouffer_score[i,j]=wilcox.test(a,b,alternative = "greater")$p.value
    }
  }
  format(round(p_val_stouffer_score, 4))

})


output$download_gtex_gt <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(gtex_d(), file)

  })
output$download_cell_gt <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(cell_grn_dt(), file)

  })


output$cell_grnplot <- downloadHandler(
  filename = function(){'cell enrich.pdf'},
  content = function(file){
    pdf(file)
    print(

      cell_grn_d()
    )
    dev.off()
  })






output$gtexplot <- downloadHandler(
  filename = function(){'gtex.pdf'},
  content = function(file){
    pdf(file)
    print(

      stouffer_gtex()
    )
    dev.off()
  })


output$hpaplot <- downloadHandler(
  filename = function(){'hpa.pdf'},
  content = function(file){
    pdf(file)
    print(

      stouffer_hpa()
    )
    dev.off()
  })
output$SCE_plot<-renderPlot({
  if(input$SCE1==0){
    return()
  }
  withProgress(message = 'Single cell Experiment Analysis...', value = 0,
               {
                  x <- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200,project = "scRNAseq")
                 
                 
                 
                 sce <<- as.SingleCellExperiment(x)
                 # is.spike <- grepl("^ERCC", rownames(sce))
                 # isSpike(sce, "ERCC") <- is.spike
                 # summary(is.spike)  
                 # is.mito <- grepl("^mt-", rownames(sce))
                 # summary(is.mito)
                 # sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
                 # head(colnames(colData(sce)))
                 # sce <- addPerCellQC(sce, subsets=list(Mito=grep("mt-", rownames(sce))))
                 # qc <- quickPerCellQC(colData(sce), sub.fields="subsets_Mito_percent")
                 # sce <- sce[,!qc$discard]
                 sce <<- computeSumFactors(sce)
                 summary(sizeFactors(sce))
                 # plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
                 #      ylab="Library size (millions)", xlab="Size factor")
                 sce <<- multiBatchNorm(sce,batch = 1)
                 dec <<- modelGeneVar(sce)
                 library(scran)
                 chosen.hvgs <<- getTopHVGs(dec, n=5000)
                 # combined <- correctExperiments(sce,PARAM=NoCorrectParam(),batch=1)
                 # 
                 library(scater)
                 set.seed(100)
                 sce <<- runPCA(sce, subset_row=chosen.hvgs)
                 sce <<- runTSNE(sce, dimred="PCA")
                 plotTSNE(sce)
               })
})



output$SCE_plot_clustering<-renderPlot({
  if(input$SCE1==0){
    return()
  }
  withProgress(message = 'Clustering Using t-SNE...', value = 0,
               {
                
                 pcs <<- reducedDim(sce, "PCA")
                 my.dist <<- dist(pcs)
                 my.tree <<- hclust(my.dist, method="ward.D2")
                 my.clusters <<- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
                 sce$cluster <<- factor(my.clusters)
                 
                 plotTSNE(sce,  colour_by="cluster",
                 ) 
                 
                 
               })
})


output$SCE_plot_silhouette <-renderPlot({
  if(input$SCE1==0){
    return()
  }
  withProgress(message =  'Calculating Silhouette Width...', value = 0,
               {
                 
                 clust.col <<- scater:::.get_palette("tableau10medium") # hidden scater colours
                 sil <<- silhouette(my.clusters, dist = my.dist)
                 sil.cols <<- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
                 sil.cols <<- sil.cols[order(-sil[,1], sil[,3])]
                 plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
                      border=sil.cols, col=sil.cols, do.col.sort=FALSE)
                 
                 
               })
})

sce_option<-reactive({
  if (input$test_use_sce=="Pairwise t-Test")
  { markers <<- findMarkers(sce, my.clusters ,test.type="t")
}
  
  else if (input$test_use_sce=="Pairwise Wilcoxon Rank Sum Test")
  { markers <<- findMarkers(sce, my.clusters ,test.type="wilcox")
}
  else (input$test_use=="Pairwise Binomial Test")
  { markers <<- findMarkers(sce, my.clusters ,test.type="binom")
  }
})
  



output$SCE_markers<-renderPrint({
  if(input$SCE1==0){
    return()
  }
  withProgress(message = 'Marker Identification...', value = 0,{
    
    # pc1 <<- reducedDim(sce, "PCA")[,1]
    # design <<- model.matrix(~pc1)
    library(limma)
    fit <<- lmFit(logcounts(sce))
    fit <<- eBayes(fit, trend=TRUE, robust=TRUE)
    topTable(fit)
    
    snn.gr <<- buildSNNGraph(sce, use.dimred="PCA")
    gr.clusters <<- igraph::cluster_fast_greedy(snn.gr)
    table(gr.clusters$membership)
    # markers <- findMarkers(sce, my.clusters, design=design, )
    sce_option()
    #pb<<-as.Seurat(sce)
    
    
   
    markers
    
    
    
    
  })
})


output$SCE_marker_display<-renderPrint({
  if(input$SCE1==0){
    return()
  }
  withProgress(message = 'Displaying the Markers...', value = 0,{
    library(limma)
    fit <<- lmFit(logcounts(sce))
    fit <<- eBayes(fit, trend=TRUE, robust=TRUE)
    topTable(fit)
    
    snn.gr <<- buildSNNGraph(sce, use.dimred="PCA")
    gr.clusters <<- igraph::cluster_fast_greedy(snn.gr)
    table(gr.clusters$membership)
    # markers <- findMarkers(sce, my.clusters, design=design, )
    sce_option()
    #pb<<-as.Seurat(sce)
 
    
    markers
    tp<-markers
    z<-c(1:length(tp))
    d = c()
    for (i in z){
      marker.set <<- markers[[i]]
      mar<-subset(marker.set,marker.set$Top<=100)
      mar_f<-subset(mar,mar$summary.logFC>=1)
      mar_f<-mar_f[,1:4]
      colnames(mar_f)[which(names(mar_f) == "p.value")] <- "p_val"
      colnames(mar_f)[which(names(mar_f) == "summary.logFC")] <- "avg_log2FC"
      mar_f$gene<-rownames(mar_f)
      mar_f$cluster<-i
      mar_f
      d<- rbind(d,mar_f)
    }
    pb.markers<<-as.data.frame(d)
  #pb.markers<<-mar_f
  d
  head(pb.markers,5)
  })
})

 
scatter_drop<-reactive({
  x <- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200,project = "scRNAseq")
  
  
  
  sce <<- as.SingleCellExperiment(x)
  
  sce <<- computeSumFactors(sce)
  summary(sizeFactors(sce))
  
  sce <<- multiBatchNorm(sce,batch = 1)
  dec <<- modelGeneVar(sce)
  library(scran)
  chosen.hvgs <<- getTopHVGs(dec, n=5000)
   
  library(scater)
  set.seed(100)
  sce <<- runPCA(sce, subset_row=chosen.hvgs)
  sce <<- runTSNE(sce, dimred="PCA")
  plotTSNE(sce)
})




output$drop_clust_scatter <- downloadHandler(
  filename = function(){'scatter_plot.pdf'},
  content = function(file){
    pdf(file)
    print( scatter_drop()
    )
    dev.off()
  })

scatter_sce<-reactive({
  
  pcs <<- reducedDim(sce, "PCA")
  my.dist <<- dist(pcs)
  my.tree <<- hclust(my.dist, method="ward.D2")
  my.clusters <<- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
  sce$cluster <<- factor(my.clusters)
  plotTSNE(sce,  colour_by="cluster",
  ) 
  
})


sce_silhouette_d<-reactive({
  clust.col <<- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <<- silhouette(my.clusters, dist = my.dist)
  sil.cols <<- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <<- sil.cols[order(-sil[,1], sil[,3])]
  plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
       border=sil.cols, col=sil.cols, do.col.sort=FALSE)
  
})

output$sce_clus_d <- downloadHandler(
  filename = function(){'scatter_plot.pdf'},
  content = function(file){
    pdf(file)
    print( scatter_sce()
    )
    dev.off()
  })

output$sce_silhouette_d <- downloadHandler(
  filename = function(){'scatter_plot.pdf'},
  content = function(file){
    pdf(file)
    print( sce_silhouette_d()
    )
    dev.off()
  })

output$feature_plot_result_sce<-renderPlot({
  if(input$SCE1==0){
    return()
  }
  scater::plotExpression(sce, features = input$feature_plot_sce 
                        )
 
})

drop_clust_table<-reactive({
  
  library(limma)
  fit <<- lmFit(logcounts(sce))
  fit <<- eBayes(fit, trend=TRUE, robust=TRUE)
  topTable(fit)
  
  snn.gr <<- buildSNNGraph(sce, use.dimred="PCA")
  gr.clusters <<- igraph::cluster_fast_greedy(snn.gr)
  table(gr.clusters$membership)
  # markers <- findMarkers(sce, my.clusters, design=design, )
  sce_option()
  #pb<<-as.Seurat(sce)
  
  
  markers
  tp<-markers
  z<-c(1:length(tp))
  d = c()
  for (i in z){
    marker.set <<- markers[[i]]
    mar<-subset(marker.set,marker.set$Top<=100)
    mar_f<-subset(mar,mar$summary.logFC>=1)
    mar_f<-mar_f[,1:4]
    colnames(mar_f)[which(names(mar_f) == "p.value")] <- "p_val"
    colnames(mar_f)[which(names(mar_f) == "summary.logFC")] <- "avg_log2FC"
    mar_f$gene<-rownames(mar_f)
    mar_f$cluster<-i
    mar_f
    d<- rbind(d,mar_f)
  }
  pb.markers<<-as.data.frame(d)
  #pb.markers<<-mar_f
  pb.markers
  
  
})


output$drop_clust_de <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(drop_clust_table(), file)
    
  })

output$sce_batch<-renderPlot({
  if(input$dropclust_batch1==0){
    return()
  }
  withProgress(message = 'Batch Correction...', value = 0,
               {
                # sce<<-SingleCellExperiment(as.matrix(data()))
                 sce1<<-batchelor::fastMNN(sce,  batch=colData(sce)[,input$meta_header_d]
                                         
                 )
                 
                 str(reducedDim(sce1, "corrected"))
                 
                 sce1 <<- runTSNE(sce1, dimred="corrected")
                 plotTSNE(sce1, colour_by="batch")
                 
                 
                 
                 
                 
                 
                 
                   })
})

batch_plot_download<-reactive({
  sce1<<-batchelor::fastMNN(sce,  batch=colData(sce)[,input$meta_header_d]
  )
  
  str(reducedDim(sce1, "corrected"))
  
  sce1 <<- runTSNE(sce1, dimred="corrected")
  plotTSNE(sce1, colour_by="batch")
  
  
})

batch_plot_download_hist<-reactive({
  
  
  cor.exp <<- assay(sce1)[input$sce_gene,]
  hist(cor.exp, xlab="Corrected expression for gene", col="#D09E2D") 
  
  
})


output$drop_clust_markers_batch<-renderPlot({
  if(input$dropclust_batch1==0){
    return()
  }
  withProgress(message = 'Corrected Expression for Gene...', value = 0,
               {
                 cor.exp <<- assay(sce1)[input$sce_gene,]
                 hist(cor.exp, xlab="Corrected expression for gene", col="#D09E2D") 
                 
})
})

sce_feature_plot<-reactive({
  scater::plotExpression(sce, features = input$feature_plot_sce 
  )
})



output$sce_feature_Plot_d <- downloadHandler(
  filename = function(){'featureplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      scater::plotExpression(sce, features = input$feature_plot_sce )
                             
      
    )
    dev.off()
  })



output$drop_clust_batch_Plot <- downloadHandler(
  filename = function(){'scatterplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      batch_plot_download()
      
      
    )
    dev.off()
  })


output$drop_clust_batch_Plot_hist <- downloadHandler(
  filename = function(){'histplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      batch_plot_download_hist()
      
      
    )
    dev.off()
  })




output$harm1<-renderPlot({
  if(input$Batch_Correction_method==0){
    return()
  }
  withProgress(message = 'Harmony Batch Correction...', value = 0,{
    
  pb <<- pb %>%
    RunHarmony(input$meta_header, plot_convergence = TRUE)
  })
})


output$harmony_plot1 <- downloadHandler(
  filename = function(){'scatterplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      pb %>%
        RunHarmony(input$meta_header, plot_convergence = TRUE)
      
      
    )
    dev.off()
  })




output$harm2<-renderPrint({
  if(input$Batch_Correction_method==0){
    return()
  }
  withProgress(message = 'Displaying Harmony Embeddings...', value = 0,{
    
    harmony_embeddings <<- Embeddings(pb, 'harmony')
    harmony_embeddings[1:5, 1:5]
  })
})
output$harm_tab <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(harmony_embeddings, file)
    
  })

output$harmony_plot2 <- downloadHandler(
  filename = function(){'scatterplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      harm3_p()
      
      
    )
    dev.off()
  })


harm3_p<-reactive({
  p1 <- DimPlot(object = pb, reduction = "harmony", pt.size = .1, group.by = input$meta_header)
  p2 <- VlnPlot(object = pb, features = "harmony_1", pt.size = .1,group.by = input$meta_header)
  plot_grid(p1,p2)
})

output$harm3<-renderPlot({
  if(input$Batch_Correction_method==0){
    return()
  }
  withProgress(message = 'Generating Feature Plot...', value = 0,{
    
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = pb, reduction = "harmony", pt.size = .1, group.by = input$meta_header)
  p2 <- VlnPlot(object = pb, features = "harmony_1", pt.size = .1, group.by = input$meta_header)
  plot_grid(p1,p2)
})
})

output$harm4<-renderPlot({
  if(input$Batch_Correction_method==0){
    return()
  }
  withProgress(message = 'Generating UMAP...', value = 0,{
    
  pb <<- pb %>%
    
    FindNeighbors(reduction = "harmony", dims = 1:20)%>%FindClusters(resolution = 0.5) %>%RunUMAP(reduction = "harmony", dims = 1:20)  %>%
     
     identity()


  umap_mat<<-pb[["umap"]]@cell.embeddings
  umap_mat<<-as.data.frame(umap_mat)
  DimPlot(pb, reduction = "umap",pt.size = 3, group.by = input$meta_header,)
  
  
  DimPlot(pb, reduction = "umap", label = TRUE, pt.size = .1)
  })
 })

harm4_p<-reactive({
  pb <<- pb %>%
    
    FindNeighbors(reduction = "harmony", dims = 1:20)%>%FindClusters(resolution = 0.5) %>%RunUMAP(reduction = "harmony", dims = 1:20)  %>%
    
    identity()
  
  
  umap_mat<<-pb[["umap"]]@cell.embeddings
  umap_mat<<-as.data.frame(umap_mat)
  DimPlot(pb, reduction = "umap",pt.size = 3, group.by = input$meta_header,)
  
  
  DimPlot(pb, reduction = "umap", label = TRUE, pt.size = .1)
})

output$harmony_plot3 <- downloadHandler(
  filename = function(){'scatterplot.pdf'},
  content = function(file){
    pdf(file)
    print(
      
      
      harm4_p()
      
      
    )
    dev.off()
  })


 metadata_add<-reactive({
  meta_data_batch<- read.csv(data1())
  pb<-AddMetaData(pb,meta_data_batch)
 })
 
 

 seurat_integrated<-reactive({
   pancreas.data <- read_rds("pancreas_expression_matrix.rds")
   metadata <- read_rds("pancreas_metadata.rds")
   pancreas <- CreateSeuratObject(pancreas.data,meta.data = metadata)
   AddMetaData(pancreas,metadata = metadata)
   
   pancreas.list <- SplitObject(pancreas, split.by = "tech")
   
   
   for (i in 1:length(pancreas.list)) {
     pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
     pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
   }
   
   reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
   pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:12)
   
   
   
   pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
   
   library(ggplot2)
   library(cowplot)
   # switch to integrated assay. The variable features of this assay are automatically
   # set during IntegrateData
   DefaultAssay(pancreas.integrated) <- "integrated"
   
   # Run the standard workflow for visualization and clustering
   pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
   pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
   pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
   p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
   #p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
               #  repel = TRUE) + NoLegend()
   #plot_grid(p1, p2)
   p1
 })
 output$meta_info2<-renderUI({
   if(input$Batch_Correction_method2==0){
     return()
   }
   pb <<- CreateSeuratObject(data(), min.cells = 3, min.features = 200,project = "scRNAseq")
  # AddMetaData(pancreas,meta.data = data1())
   
   meta_d<<-data1()
   pb<<-AddMetaData(pb,meta_d)
   choices_meta<<-colnames(pb@meta.data)
   selectInput("meta_header","Select Parameter for Batch Correction",choices = choices_meta)
 })
 output$seurat_integrated_plot<-renderPlot({
   if(input$seurat_int1==0){
     return()
   }
   withProgress(message = 'Batch Correction...', value = 0,{
   # pancreas.data <- read_rds("pancreas_expression_matrix.rds")
   # metadata <- read_rds("pancreas_metadata.rds")
   # pb <<- CreateSeuratObject(data(), min.cells = 3, min.features = 200,project = "scRNAseq")
   # # AddMetaData(pancreas,meta.data = data1())
   # 
   # meta_d<<-data1()
   # pb<<-AddMetaData(pb,meta_d)
   #pb@meta.data<<-meta_d
   pb.list <<- SplitObject(pb, split.by = input$meta_header)
   
   
   for (i in 1:length(pb.list)) {
     pb.list[[i]] <- NormalizeData(pb.list[[i]], verbose = FALSE)
     pb.list[[i]] <- FindVariableFeatures(pb.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
   }
   
   #reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
   pb.anchors <<- FindIntegrationAnchors(object.list = pb.list, dims = 1:30)
   
   
   
   pb.integrated <<- IntegrateData(anchorset = pb.anchors, dims = 1:30)
   
   library(ggplot2)
   library(cowplot)
   # switch to integrated assay. The variable features of this assay are automatically
   # set during IntegrateData
   DefaultAssay(pb.integrated) <<- "integrated"
   
   # Run the standard workflow for visualization and clustering
   pb<<- ScaleData(pb.integrated, verbose = FALSE)
   pb<<- RunPCA(pb, npcs = 30, verbose = FALSE)
   pb <<- RunUMAP(pb, reduction = "pca", dims = 1:30)
   p1 <- DimPlot(pb, reduction = "umap", group.by = input$meta_header)
   #p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
   #  repel = TRUE) + NoLegend()
   #plot_grid(p1, p2)
   p1
   })
 })
 


})
