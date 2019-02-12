require(broom)

cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")

### Slignshot ##

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


runSDS  <- function(object,reduction.type='dm',dim.use=1:2,groups=NULL, start.clus=NULL,end.clus=NULL){
  rd <-GetCellEmbeddings(object,reduction.type,dim.use)
  cl <- object@ident
  object@misc[['sds']] <-  list("dr"=reduction.type,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus))
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


plotPseudoTime = function(object,groupby){

  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
                                                d <- as.data.frame(c$s[c$ord,seq_len(2)])
                                                d$curve<-x
                                                return(d)
                                                })
                                               ) 

  p=FetchData(object,append(c('DM1','DM2'),groupby)) %>%
    ggplot(.,aes(x=DM1,y=DM2))+geom_point(aes(color=!!sym(groupby))) + theme(legend.position="top") + guides(col = guide_legend(nrow = 2))+
    geom_path(aes(DM1, DM2,linetype=curve),curved,size=1)
p
  
  
    
  }
  
  
##




