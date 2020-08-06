#' @title The default pipline for processing scRNA-seq by Seurat
#'
#' @description Analysis the single cell data according to the default parameters of Seurat. Users can adjust the parameters according to their needs.
#'
#' @description This step can be skip, if user can provided the processed data.
#' @description The processed data should finish several steps by using Seurat: Normalized data, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors, FindClusters and RunUMAP.
#'
#'
#' @param data.dir: Path to the matrix from CellRanger reults, such as "~/sample_CellRanger_result/outs/filtered_feature_bc_matrix"
#' @param PC: Number of PC used for demension reduction and clustering. default:40
#' @param resolution: Choose appropriate resolution value to adjust the number of clusters. default:0.6
#' @param mt.cut_off: Filter out the cells with a high percentage of MT-genes. default:20(%) (Only apply to samples for human)
#' @param min_nFeature.cut_off: Filter low quality cell by number of Features in each cells. default:200
#' @param sample_name: Samples name of scRNA-seq data. default:"scRNA-seq"
#' @param data_type: Type of scRNA-seq data. Temporary, only "Expression_matrix" and "TotalSeq" data are supported
#' @param filter_doublet: Whether use DoubletFinder to remove Doublet in scRNA-seq; Attention, DoubletFinder performed not pretty well in detect the doublet of same type of cells.
#' @export
#'
#'
Seurat_Pipline<-function(data.dir,
                         PC=40,
                         resolution=0.6,
                         mt.cut_off=20,
                         min_nFeature.cut_off=200,
                         sample_name="scRNA-seq",
                         data_type="Expression_matrix",
                         filter_doublet=TRUE
){
  #Load dependent packages
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(glue)

  #Load data


  object.data<-Read10X(data.dir)


  if( data_type=="Expression_matrix" ){
    object <- CreateSeuratObject(counts = object.data, project = sample_name, min.cells = 3, min.features = 200)

  } else if( data_type=="TotalSeq" ){
    object <- CreateSeuratObject(counts = object.data[[1]], project = sample_name, min.cells = 3, min.features = 200)
  } else{
    print("Error:Please choose from Expression_matrix or TotalSeq ")
  }

  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

  # QC
  if(min_nFeature.cut_off==0){
    if(median(object$nFeature_RNA)<1000){
      object <- subset(object, subset = nFeature_RNA > 200  & percent.mt < mt.cut_off)

    } else{
      object <- subset(object, subset = nFeature_RNA > 500  & percent.mt < mt.cut_off)
    }
  }else {
    object <- subset(object, subset = nFeature_RNA > min_nFeature.cut_off  & percent.mt < mt.cut_off)
  }


  # Normalize
  object<- NormalizeData(object,normalization.method = "LogNormalize", scale.factor = 10000)
  # Feature Selection
  object <- FindVariableFeatures(object,
                                 selection.method = "vst",
                                 nfeatures = 2000)
  object <- ScaleData(object,features=VariableFeatures(object))
  # demension reduction
  object <- RunPCA(object,npcs =PC,verbose = F,features=VariableFeatures(object))

  # Clustering
  object <- FindNeighbors(object, dims =1:PC)
  object <- FindClusters(object, resolution = resolution)
  # Visualization
  object <- RunUMAP(object, dims = 1:PC)

  # Rename the ID start from 1
  new_idents<-as.numeric(levels(object))+1
  names(new_idents)<-levels(object)
  object<-RenameIdents(object,new_idents)
  object$First_time_unsupervised_clustering<-Idents(object)
  if(filter_doublet==TRUE){
    object<-run_DoubletFinder(object,PC=PC,resolution = resolution ,sample_name =sample_name ,removeDoublet = T)
  }

  return(object)

}



run_DoubletFinder<-function(object,PC=40,resolution=0.6,sample_name="scRNA-seq",removeDoublet=TRUE){
  require(DoubletFinder) # load DoubletFinder
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(glue)
  sweep.data <- paramSweep_v3(object,PCs=1:PC)
  sweep.stats <- summarizeSweep(sweep.data, GT = FALSE)
  bcmvn= find.pK(sweep.stats) # Find the best pK value
  homotypic.prop=modelHomotypic(object@meta.data[,paste0("RNA_snn_res.",resolution)])
  nExp_poi=round(0.075*length(object$orig.ident))
  nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
  object=doubletFinder_v3(object, PCs = 1:PC, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),
                          nExp = nExp_poi.adj, reuse.pANN = FALSE) # Identify Doublet

  object@meta.data$DF_hi.lo<- object@meta.data[,dim(object@meta.data)[2]]
  Idents(object ) <- "DF_hi.lo"
  if(removeDoublet==TRUE){
    object<-subset(x = object, idents="Singlet") # Remove Doublet

    # Repeat these steps
    # Normalize
    object<- NormalizeData(object,normalization.method = "LogNormalize", scale.factor = 10000)
    # Feature Selection
    object <- FindVariableFeatures(object,
                                   selection.method = "vst",
                                   nfeatures = 2000)
    object <- ScaleData(object,features=VariableFeatures(object))
    # demension reduction
    object <- RunPCA(object,npcs =max(PC),verbose = F,features=VariableFeatures(object))

    # Clustering
    object <- FindNeighbors(object, dims =1:PC)
    object <- FindClusters(object, resolution = resolution)
    # Visualization
    object <- RunUMAP(object, dims = 1:PC)

    # Rename the ID start from 1
    new_idents<-as.numeric(levels(object))+1
    names(new_idents)<-levels(object)
    object<-RenameIdents(object,new_idents)

  }
  return(object)
}




