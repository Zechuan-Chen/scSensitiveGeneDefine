# scSensitiveGeneDefine

---
## **Attention:** 
All the code has been integrated into the R packages "scSensitiveGeneDefine";

This repository will be renamed as "scSensitiveGeneDefine"

## Description:
`scSensitiveGeneDefine` is a R package that define the sensitive genes in single-cell RNA sequencing data.

`scSensitiveGeneDefine` is build based on Seurat(>= 3.0.1)(https://satijalab.org/seurat/); DoubletFinder(>= 2.0.3)(https://github.com/chris-mcginnis-ucsf/DoubletFinder); entropy>=1.2.1; ;dplyr(>=1.0.0); All of these four dependent packages are R package.

`scSensitiveGeneDefine` intend to publish on BMC Bioinformatics.

## Installation(in R/Rstudio)
devtools::install_github("Zechuan-Chen/scSensitiveGeneDefine")

## Dependencies
`scSensitiveGeneDefine` requires the following R packages:

 - Seurat (>=3.0.1)
 - DoubletFinder (>=2.0.3)
 - entropy (>=1.2.1)
 - dplyr (>=1.0.0)
 - NOTE:The version of these depend packages are temporary.

Example code for `scSensitiveGeneDefine`

```
object<-runSeurat(data.dir="~/outs/filtered_feature_bc_matrix/",
                       sample_name = "scRNA-seq Sample 1",
                       PC = 40,
                       resolution = 0.6,
                       mt.cut_off = 20,
                       min_nFeature.cut_off = 200,
                       data_type = "Expression_matrix",
                       filter_doublet = T,
					   algorithm=1)
                       
# The processed object also can be provided by user!

HVG_Anno<-HVG_Statistic(object)
SensitiveGene<-GetSensitivegene(object,min_nClusters = "Default",HVG_Anno = HVG_Anno)
object<-ReSelectVariableFeatures(object,SensitiveGene = SensitiveGene)
object<-ReClustering(object,PC = 40,resolution = 0.6,algorithm=1)

# Evaluate the clustering result (If you have the grount-truth labels)

ECA_value<-ECA(object,Ground_truth_label = label1,Generated_label = label2)
ECP_value<-ECP(object,Ground_truth_label = label1,Generated_label = label2)
```




