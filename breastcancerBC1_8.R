# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 3 - EVA breast tumor vs normal
# code developed by Emily Davis

sink = (file = "rout.log")
library(monocle)
library(GSReg)
print("libraries loaded")

###################################################
### code chunk 1: load biscuit imputed data BC 1-8
###################################################

load("breast_cancer_imputed_cdsBC1_8.rda")
print("cds loaded")
load("hallmark_pathways.rda")
print("pathways loaded")

###################################################
### code chunk 4: eva
###################################################

# summary 
table(pData(cds)$Celltype, pData(cds)$tissue)

do_pathways <- function(pathway, type) {
  
  # subset tumor and normal tissues
  test <- cds[,(pData(cds)$tissue=="TUMOR" | pData(cds)$tissue=="NORMAL")]
  testc <- test[,grep(type, pData(test)$Celltype) ]
  pData(testc)$compare <- ifelse(pData(testc)$tissue == "NORMAL", 0,1)
  #celltype1 =0; celltype2 =1 in pData(testc)$compare
  # get expression matrix
  exprsdata <- as.matrix(exprs(testc))
  rownames(exprsdata) <- fData(testc)$gene_short_name
  #classic eva
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway, phenotypes=as.factor(pData(testc)$compare)) 
  
  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)
  
  dd$metric <- "kendall_tau"
  dd$celltype <- type
  dd$group1 <- "Normal"
  dd$group2 <- "Tumor"
  return(dd)
}

a <- do_pathways(pathway = hallmark_pathways, type = "T:CD8")
print("comparison 1 done")
b <- do_pathways(pathway = hallmark_pathways, type = "MONOCYTE:")
print("comparison 2 done")
c <- do_pathways(pathway = hallmark_pathways, type = "NK:")
print("comparison 3 done")
d <- do_pathways(pathway = hallmark_pathways, type = "T:CD4")


