#' @title Reclustering after removing sensitive genes
#'
#' @description After removing the sensitive gene, reclusering the single-cell data.
#'
#' @param object: Seurat object, with HVG_Statistic, GetSensitiveGene and ReSelectVariableFeatures finished
#' @param PC: Number of PC used for demension reduction and clustering. default:40
#' @param resolution: Choose appropriate resolution value to adjust the number of clusters. default:0.6
#' @param algorithm: Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @export
#'

ReClustering<-function(object,
                       PC=40,
                       resolution=0.6,
                       algorithm=1
){
 
  object <- ScaleData(object,features=VariableFeatures(object))
  # demension reduction
  object <- RunPCA(object,npcs =max(PC),verbose = F,features=VariableFeatures(object))

  # Clustering
  object <- FindNeighbors(object, dims =1:PC)
  object <- FindClusters(object, resolution = resolution,algorithm=algorithm)
  # Visualization
  object <- RunUMAP(object, dims = 1:PC)

  # Rename the ID start from 1
  new_idents<-as.numeric(levels(object))+1
  names(new_idents)<-levels(object)
  object<-RenameIdents(object,new_idents)

  object$Reclustering<-Idents(object)
  return(object)
}

