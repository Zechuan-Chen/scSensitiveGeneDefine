#' @title Screening for sensitive genes
#' @description Screening for sensitive genes based on the CV Ranking and Shannon Index
#'
#'
#' @param object: The seurat object, with HVG_Statistic finished.
#' @param min_nClusters: nClusters of sensitive genes should be detect as HVGs at least,default:"Default"(Half of nCluster)
#' @param min_nCell: nCells of sensitive genes should be detect as HVGs as least,default:0
#' @param HVG_Anno: Statistics result on the high variable genes in each subclusterd
#' @export
GetSensitivegene<-function(object,
                           min_nClusters="Default",
                           min_nCell=0,
                           HVG_Anno
){
  require(entropy)
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(glue)
  if(min_nClusters=="Default"){
    min_nClusters<-floor(length(levels(object))/2)
  }

  pre_Senstive_gene.statistic<-HVG_Anno$HVG_Type.statistic[HVG_Anno$HVG_Type.statistic$nCluster > min_nClusters & HVG_Anno$HVG_Type.statistic$nCell > min_nCell,]
  pre_Senstive_gene<-rownames(pre_Senstive_gene.statistic)
  avg_expression<-AverageExpression(object ,assays = "RNA",features = pre_Senstive_gene)
  data<-apply(avg_expression$RNA,1,function(x) entropy(x) )
  data<-as.data.frame(data)
  data$gene<-rownames(data)
  data<-data[order(data$data),]
  emtropu_value<-data
  emtropu_value$nCluster<-pre_Senstive_gene.statistic[rownames(data),]$nCluster
  emtropu_value$nCell<-pre_Senstive_gene.statistic[rownames(data),]$nCell
  colnames(emtropu_value)<-c("entropy_value","gene","nCluster","nCell")
  SensitiveGene<-data[data$data> median(data$data),]
  colnames(SensitiveGene)<-c("entropy_value","gene")
  return(SensitiveGene)
}
