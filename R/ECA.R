#' @title Caculate the ECA Value
#'
#' @description Evaluate the clustering result by useing ECA(entropy of accuracy)
#'
#' @param object: Seurat object, with ground-truth labels and clustering labels
#' @param Ground_truth_label: The real cell type labels, saved in meta.data of Seurat object
#' @param Generated_label: The clustering labels or others labels, saved in meta.data of Seurat object
#' @export
ECA<-function(object,
              Ground_truth_label,
              Generated_label
){
  Numerator<-c()
  x<-object@meta.data
  Denominator<- x[,Generated_label]%>% table() %>% names() %>% length()
  rownames(x)<-c(1:dim(x)[1])
  x[,Ground_truth_label]<-as.vector(x[,Ground_truth_label])
  x[,Generated_label]<-as.vector(x[,Generated_label])
  for(i in names(table(x[,Generated_label]))){
    level<-x[x[,Generated_label]==i,][,Ground_truth_label] %>% table() %>% names()
    xi<-x[x[,Generated_label]==i,]
    for(j in level){
      xij<-xi[xi[,Ground_truth_label]==j,]
      p_xij<-dim(xij)[1]/dim(xi)[1]
      tmp_Numerator<-p_xij*log(p_xij)
      Numerator<-c(Numerator,tmp_Numerator)
    }
  }
  epa= -sum(Numerator)/Denominator
  return(epa)
}
