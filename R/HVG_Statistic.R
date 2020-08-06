#' @title Statistics on high variable genes of each subclusters
#' @description Statistics on the high variable genes in each subclusterd, divided by the First-time unsupervised clustering
#' @description Only need to provide the label of First-time unsupervised clustering.
#'
#' @param object: The seurat object, produced by Seurat_Pipline Function or provided by User
#' @param First_time_unsupervised_clustering_label: The labels of Unsupervised clustering results
#' @export
HVG_Statistic<-function(object,First_time_unsupervised_clustering_label="First_time_unsupervised_clustering"){
  require("entropy")
  Idents(object)<-First_time_unsupervised_clustering_label
  label1<-levels(object) # Unsupervised clustering labels
  subtype_cells<-list() # Each Cluster
  Variable_list<-list() # Variable gene list
  Common_HVG<-c() # number of HVGs in each cluster
  meta.features<-list() # meta data of features selection
  for(i  in 1:length(label1)){
    subtype_cells[[label1[i]]] <- subset(object,idents=label1[i])
    subtype_cells[[label1[i]]] <- FindVariableFeatures(subtype_cells[[label1[i]]], selection.method = "vst", nfeatures = 2000)
    Variable_list[[label1[i]]] <- VariableFeatures(subtype_cells[[label1[i]]])
    meta.features[[label1[i]]] <- subtype_cells[[label1[i]]]@assays$RNA@meta.features
    if(i != 1){
      Common_HVG<- intersect(Common_HVG,Variable_list[[label1[i]]])
    }else{
      Common_HVG=Variable_list[[label1[i]]]
    }
  }
  Common_HVG.statistic<-data.frame(row.names=label1,nHVG=rep(NA,length(label1)))
  emtropu_value<-c()
  for(j  in 1:length(label1)){
    emtropu_value<-c(emtropu_value,length(Variable_list[[label1[j]]]))
  }

  Common_HVG.statistic$nHVG<-emtropu_value
  Total_HVG<-c()
  for(k  in 1:length(label1)){
    Total_HVG<-c(Total_HVG,Variable_list[[label1[k]]])
  }
  HVG_Type<-unique(Total_HVG)
  HVG_Type.statistic<-data.frame(row.names=HVG_Type,nCluster=rep(NA,length(HVG_Type)),nCell=rep(NA,length(HVG_Type)))
  x<-as.data.frame(table(Total_HVG))
  rownames(x)<-x$Total_HVG
  HVG_Type.statistic$nCluster<-x[HVG_Type,]$Freq


  for(l in 1:length(HVG_Type)){
    nCell=0
    rm(i)
    for(i  in 1:length(label1)){
      if(HVG_Type[l] %in% Variable_list[[label1[i]]]){
        nCell=nCell+dim(subtype_cells[[label1[i]]]@meta.data)[1]
      }
    }
    HVG_Type.statistic[l,2]<-nCell
  }
  HVG_Anno<-list(Variable_list=Variable_list,
                 Common_HVG.statistic=Common_HVG.statistic,
                 HVG_Type.statistic=HVG_Type.statistic,
                 meta.features=meta.features)
  return(HVG_Anno)
}
