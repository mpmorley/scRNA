require(broom)
library(slingshot)

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


RunDiffusion <- function(
  object,
  dims = 1:5,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  max.dim = 2L,
  q.use = 0.01,
  reduction.name = "dm",
  reduction.key = "DM_",
  ...
) {
  if (!is.null(x = dims) || is.null(x = features)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
  } else {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  }
  
  data.dist <- dist(data.use)
  data.diffusion <- data.frame(destiny::DiffusionMap(data = as.matrix(data.dist),n_eigs = max.dim)@eigenvectors)
  
  colnames(x = data.diffusion) <- paste0(reduction.key, 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <-  rownames(data.use)
 # for (i in 1:max.dim) {
#    x <- data.diffusion[, i]
 #   x <- MinMax(data = x, min = quantile(x = x, probs = q.use), 
  #              quantile(x = x, probs = 1 - q.use))
  #  data.diffusion[, i] <- x
  #}
  
  assay <- DefaultAssay(object = object[[reduction]])
  
  dm.reduction <- CreateDimReducObject(
    embeddings = as.matrix(data.diffusion),
    key = reduction.key,
    assay = assay
  )
  
  
  object[[reduction.name]] <- dm.reduction
  
  
  
  #object <- LogSeuratCommand(object = object)
  return(object)
}


################################
#   Slingshot 
#
###################################
plotCurveHeatmaps <- function(object=NULL,curve=NULL,filename='heatmap.png',n=25){
    c = sym(curve)
    cells = object@meta.data %>% tibble::rownames_to_column('cellid') %>% arrange(desc(!!c)) %>% filter(!is.na(!!c))
    genes = object@misc$sds$dge[[curve]] %>% arrange(p.value) %>% head(n) %>% pull(gene)
    FetchData(object,genes,cells$cellid,use.scaled=T) %>% t(.) %>%
      NMF::aheatmap(.,Colv=NA,distfun='pearson',scale='row',annCol=cells$var_celltype,annColors = list(X1=cpallette), 
                    filename=filename)

}


plotCurveDGEgenes <- function(object=NULL,curve=NULL,n=25,reduction.use='dm'){
  genes = object@misc$sds$dge[[curve]] %>% arrange(p.value) %>% head(n) %>% pull(gene)
  plot_grid(  plotlist = FeaturePlot(scrna.sub,genes,reduction.use = reduction.use,cols.use = c('grey','purple'),do.return = T))
  
}




runSDS  <- function(object,reduction='dm',groups=NULL, start.clus=NULL,end.clus=NULL){
  rd <- Embeddings(object,reduction)
  cl <- Idents(object = object)
  object@misc[['sds']] <-  list("dr"=reduction,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus))
  ps <- slingPseudotime(object@misc[['sds']]$data)
  object@meta.data[,colnames(ps)] <- ps 
  return(object)
}


runSDSDGE <- function(object){
  DGE <- list()
  for(c in names(object@misc$sds$data@curves)){
    object@misc$sds$dge[[c]] <- FetchData(object,append(object@var.genes, c,0),use.scaled = T ) %>% tidyr::gather(gene,signal, -one_of(c)) %>% dplyr::rename(curve = 1) %>%
      tidyr::nest(-gene) %>% 
     mutate(
       fit = purrr::map(data, ~ gam(signal ~ lo(curve), data = .x)),
        tidied = purrr::map(fit, tidy)
      ) %>% 
      tidyr::unnest(tidied) %>% 
      filter(term !='Residuals')
  }
  return(object)

}


plotPseudoTime = function(object,groupby,reduction.use='DM'){

  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
                                                d <- as.data.frame(c$s[c$ord,seq_len(2)])
                                                d$curve<-x
                                                return(d)
                                                })
                                               ) 

  
  data <- Embeddings(object = object[[reduction]])
  
  p=FetchData(object,groupby) %>% 
    ggplot(.,aes(x=DM_1,y=DM_2))+geom_point(aes(color=!!sym(groupby))) + theme(legend.position="top") + guides(col = guide_legend(nrow = 2))+
    geom_path(aes(DM_1, DM_2,linetype=curve),curved,size=1)
p
  
  
    
  }
  
  
  
  




