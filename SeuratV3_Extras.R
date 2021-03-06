require(broom)
require(plotly)
#require(slingshot)
require(dplyr)
devtools::install_github("Morriseylab/ligrec")



cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")


processExper <- function(dir,name,org='mouse',files,ccscale=F,filter = T){
  try(if(length(files)==0) stop("No files"))
 
  
  if(length(files)==1){
    # Load the dataset
    inputdata <- Read10X(data.dir =files[1])
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name)
    # Initialize the Seurat object with the raw (non-normalized data).  
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200,project = name)
  }else{
    #Initialize the first object with the raw (non-normalized data) and add rest of the data 
    inputdata <- Read10X(data.dir =files[1])
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name, '-rep1')
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200, project = name)
    #cat('Rep1', length(object@cell.names), "\n")
    for(i in 2:length(files)){
      tmp.data <- Read10X(data.dir =files[i])
      colnames(tmp.data) <- paste0(colnames(tmp.data), '-',name, '-rep',i)
      
      tmp.object <- CreateSeuratObject(counts= tmp.data, min.cells = 10, min.features = 200, project = name)
     # cat('Rep', i, ": ", length(tmp.object@cell.names), "\n", sep="")
      object <- merge(object, tmp.object, do.normalize = FALSE, min.cells = 0, min.features = 0)
    }
    #cat("merged: ", length(object@cell.names), "\n", sep="")
  }
  

  
  if(org=='mouse'){
    mito.features <- grep(pattern = "^mt-", x = rownames(x = object), value = TRUE)
  }else{
    mito.features <- grep(pattern = "^MT-", x = rownames(x = object), value = TRUE)
  }
  
  
  percent.mito <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts")[mito.features,
                                                                                     ])/Matrix::colSums(x = GetAssayData(object = object, slot = "counts"))
  
  # The [[ operator can add columns to object metadata, and is a great place
  # to stash QC stats
  object[["percent.mito"]] <- percent.mito
  VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
          ncol = 3)
  
  if(filter){
    #Using a median + 3 MAD cutoff for high genes. 
    object <- subset(object, subset = nFeature_RNA > 200 & percent.mito < 0.05 & nFeature_RNA < median(object$nFeature_RNA) + 3*mad(object$nFeature_RNA) )
    
  }
  
  #normalize data
  object <- NormalizeData(object = object)
  
  #detection of variable genes
  #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
  object <-FindVariableFeatures(object = object, 
                                selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  
  if(ccscale==T){
    if(org=='human'){
      #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
      object <- CellCycleScoring(object = object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }else{
      m2h <- readr::read_csv(mouseorthologfile)
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

ClusterDR <-function(object,npcs=50, maxdim='auto',k=30){
  object <- RunPCA(object = object, npcs = npcs, verbose = FALSE)
  
  if(maxdim=='auto'){
    object <- JackStraw(object = object, num.replicate = 100,dims = npcs)
    object <- ScoreJackStraw(object = object,dims=1:npcs)
    dim <- object@reductions$pca@jackstraw$overall.p.values %>% 
      as.data.frame(.) %>% 
      mutate(adj = p.adjust(Score,method='bonferroni')) %>% 
      filter(adj <0.05) %>% 
      summarise(max=max(PC)) %>% 
      pull(max)
    
  } else {
    dim<-maxdim
  }
  print(dim)
  object <- RunTSNE(object = object, reduction = "pca",dims = 1:dim)
  object <- RunUMAP(object = object, reduction = "pca", n.neighbors = k,n.components = 3,dims = 1:dim)
  object <- RunDiffusion(object = object,dims=1:dim)
  object <- FindNeighbors(object = object,dims=1:dim,k.param = k)
  object <- FindClusters(object = object)
  object$var_cluster <- object@active.ident
  object@misc[["findallmarkers"]] <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  object
}


getMaxDim <- function(object){
  
  object@reductions$pca@jackstraw$overall.p.values %>% 
    as.data.frame(.) %>% 
    mutate(adj = p.adjust(Score,method='bonferroni')) %>% 
    filter(adj <0.05) %>% 
    summarise(max=max(PC)) %>% 
    pull(max)
  
  
}

getClusterMarkers <- function(object,cluster=0){
  
  object@misc[['findallmarkers']] %>% filter(cluster==!!cluster)
  
}



HeatMapTopGenes <- function(object,nfeatures=10){
  if(nfeatures > 50){
    print('Too many features, choose < 50')
    return()
  }
 
  object@misc[['findallmarkers']] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene) %>%
  DoHeatmap(object = object, features = .) + NoLegend()
  
}

DotPlotTopGenes <- function(object,nfeatures=3){
  if(nfeatures > 50){
    print('Too many features, choose < 10')
    return()
  }
  
  object@misc[['findallmarkers']] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene) %>%
    DotPlot(object = object, features = .)
  
}


RunDiffusion <- function(
  object,
  dims = 1:5,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  max.dim = 3L,
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
  
  data.dist <- parallelDist::parDist(data.use)
  data.diffusion <- data.frame(destiny::DiffusionMap(data = as.matrix(data.dist)+1,n_eigs = max.dim)@eigenvectors)
  
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

make3dPlot <- function(object,groupby,reduction='dm',colors=NULL){
  dims=1:3
  dims <- paste0(Key(object = object[[reduction]]), dims)
  data <- FetchData(object = object, vars = c(dims,groupby))
  
  if(is.factor(data[,groupby])){
    colors=cpallette
  }
  plot_ly(data, x=~get(dims[1]), y=~get(dims[2]), z=~get(dims[3]),colors=colors,color=~get(groupby),size=.5 ) %>%
    add_markers()  
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



runSDS <- function(object,reduction='dm',groups=NULL, start.clus=NULL,end.clus=NULL){
 print('Please use the runSlingshot Function') 
}

runSlingshot  <- function(object,reduction='dm',groups=NULL, start.clus=NULL,end.clus=NULL){
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





plotPseudoTime = function(object,groupby,reduction='dm',dims=1:2){

  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
                                                d <- as.data.frame(c$s[c$ord,seq_len(2)])
                                                d$curve<-x
                                                return(d)
                                                })
                                               ) 

  dims <- paste0(Key(object = object[[reduction]]), dims)

  p=FetchData(object = object, vars = c(dims,groupby)) %>%
    ggplot(.,aes_string(x=dims[1],y=dims[2]))+geom_point(aes(color=!!sym(groupby))) + 
    theme(legend.position="top") + 
    guides(col = guide_legend(nrow = 2)) +
    geom_path(aes_string(dims[1], dims[2],linetype="curve"),curved,size=1)
p
  
  
    
  }
#####################################################################################
#
#
#####################################################################################

ligrec <- function(object,grp.var='ident',org,perc=30){
  #get grouping variable
  var=as.character(grp.var)
  #genes=fread("data/ligrecgenes.txt",header = TRUE)       
  if(org=="mouse"){
    data('Mm_PairsLigRec',package="ligrec")
    rl = mm
  }else if(org=="human"){
    data('Hs_PairsLigRec',package="ligrec")
    rl=hs
    }

  genes <- intersect(rownames(GetAssayData(object = object, slot = "counts",assay='RNA')), unique(c(as.character(rl$ligand),as.character(rl$receptor))))
  
  #For all unique genes in the ligrec list, get their expression value for all cells and the groups the cells belong to
  my.data <- cbind(FetchData(object,c(var)), FetchData(GetAssay(object = object,assay='RNA'),genes,slot="counts"))
  colnames(my.data)[1]= "clust"
  perc=perc/100
  result=data.frame()
  res=data.frame()
  #loop over each cluster to find pairs
  for(i in 1:(length(levels(my.data$clust)))){
    for(j in 1:(length(levels(my.data$clust)))){
      #from the large martix, subselect receptor and lig subgoups (if i=1 and j=2, keep cells in grps 1 and 2)
      test=my.data[my.data$clust==levels(my.data$clust)[i] | my.data$clust==levels(my.data$clust)[j],]
      #Subselect genes in receptor list in cells in rec subgroup (say 1)
      R_c1=test[test$clust==levels(my.data$clust)[i] ,(colnames(test) %in% rl$receptor)]
      #Subselect genes in ligand list in cells in lig subgroup (say 2)
      L_c2=test[test$clust==levels(my.data$clust)[j] , (colnames(test) %in% rl$ligand)]
      if(nrow(R_c1)!=0 &nrow(L_c2)!=0){
        #keep genes that are expressed in more than user-input percent of the cells
        keep1 = colSums(R_c1>0)>=perc*dim(R_c1)[1]
        keep2 = colSums(L_c2>0)>=perc*dim(L_c2)[1]
        R_c1=R_c1[,keep1]
        L_c2=L_c2[,keep2]
        #get list of lig-rec pairs
        res=rl[(rl$ligand %in% colnames(L_c2)) & (rl$receptor %in% colnames(R_c1)),]
      }else{}
      if(nrow(res)!=0){
        res$Receptor_cluster=levels(my.data$clust)[i]
        res$Lig_cluster=levels(my.data$clust)[j]
        result=rbind(result,res)
      }else{result=result}
    }
  }
  # get final list of all lig-rec pairs
  #result=result[result$Receptor_cluster!=result$Lig_cluster,]
  
  object@misc[['ligrec']] <- result
  object
}

  




