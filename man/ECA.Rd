% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECA.R
\name{ECA}
\alias{ECA}
\title{Caculate the ECA Value}
\usage{
ECA(object, Ground_truth_label, Generated_label)
}
\arguments{
\item{object:}{Seurat object, with ground-truth labels and clustering labels}

\item{Ground_truth_label:}{The real cell type labels, saved in meta.data of Seurat object}

\item{Generated_label:}{The clustering labels or others labels, saved in meta.data of Seurat object}
}
\description{
Evaluate the clustering result by useing ECA(entropy of accuracy)
}
