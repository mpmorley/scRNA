
processExper <- function(dir,name,org='mouse',files,ccscale=F){
  try(if(length(files)==0) stop("No files"))
  
    # Load the dataset
    inputdata <- Read10X(data.dir =files[1])
    colnames( inputdata) <- paste0(colnames(inputdata), '_',name)
    
    # Initialize the Seurat object with the raw (non-normalized data).  
    object <- CreateSeuratObject(counts = inputdata, min.cells = 10, min.features = 200,project = name)
    
    
    if(org=='mouse'){
      mito.features <- grep(pattern = "^mt-", x = rownames(x = object), value = TRUE)
    }else{
      mito.features <- grep(pattern = "^Mt-", x = rownames(x = object), value = TRUE)
    }
    
    
    percent.mito <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts")[mito.features,
                                                                                     ])/Matrix::colSums(x = GetAssayData(object = object, slot = "counts"))
    
    # The [[ operator can add columns to object metadata, and is a great place
    # to stash QC stats
    object[["percent.mito"]] <- percent.mito
    VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
            ncol = 3)
    

  #normalize data
    object <- NormalizeData(object = object)
  
  #detection of variable genes
  #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
    object <-FindVariableFeatures(object = object, 
                               selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
    
    if(ccscale==T){
      if(org=='human'){
        #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
        object <- CellCycleScoring(object = object, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
      }else{
        m2h <- read_csv(mouseorthologfile)
        cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
        cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
        object <- CellCycleScoring(object = object, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
      }
      #Scaling the data and removing unwanted sources of variation
      object <- ScaleData(object = object, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
    }else{
      object <- ScaleData(object = object, vars.to.regress = c("nUMI", "percent.mito"))
    }
  return(object)
    
    }

