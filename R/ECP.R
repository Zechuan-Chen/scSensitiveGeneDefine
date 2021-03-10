#' @title Caculate the ECP Value
#'
#' @description Evaluate the clustering result by useing ECP(entropy of purity)
#'
#' @param object: Seurat object, with ground-truth labels and clustering labels
#' @param Ground_truth_label: The real cell type labels, saved in meta.data of Seurat object
#' @param Generated_label: The clustering labels or others labels, saved in meta.data of Seurat object
#' @export
ECP<-function(x=object,
Ground_truth_label,
Generated_label
){
  Numerator<-c()
  Denominator<- x[,Ground_truth_label]%>% table() %>% names() %>% length()
  rownames(x)<-c(1:dim(x)[1])
  x[,Ground_truth_label]<-as.vector(x[,Ground_truth_label])
  x[,Generated_label]<-as.vector(x[,Generated_label])
  for(i in names(table(x[,Ground_truth_label]))){
    level<-x[x[,Ground_truth_label]==i,][,Generated_label] %>% table() %>% names()

    for(j in level){
      xi<-x[x[,Generated_label]==j,]
      xij<-xi[xi[,Ground_truth_label]==i,]
      p_xij<-dim(xij)[1]/dim(xi)[1]
      tmp_Numerator<-p_xij*log(p_xij)
      Numerator<-c(Numerator,tmp_Numerator)
    }
  }
  eca= -sum(Numerator)/Denominator
  return(eca)
}
