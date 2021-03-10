#' @title Remove sensitive genes from variable Features
#'
#' @description Only removing the sensitive genes from the variable features, not transcriptome. In this way, the unsupervised clustering result will be improved with transcript no change.
#' @param object: Seurat object, with HVG_Statistic and GetSensitiveGene finished
#' @param SensitiveGene: SensitiveGene List of the scRNA-seq data. The sensitive gene can be get by running "GetSensitiveGene"
#' 
#' 
#' @export
#' 
ReSelectVariableFeatures<-function(object,SensitiveGene){
  matrix<-object@assays$RNA@meta.features
  matrix<-matrix[order(matrix$vst.variance.standardized,decreasing=T),]
  matrix<-matrix[setdiff(rownames(matrix),SensitiveGene$gene),]
  new_Features<-rownames(matrix)[1:2000]
  object@assays$RNA@meta.features$vst.variable<-FALSE
  object@assays$RNA@meta.features[new_Features,]$vst.variable<-TRUE
  VariableFeatures(object)<-new_Features
  return(object)
}

