
###################################################
### code chunk: functions
###################################################

splat_eva <- function(mat, pheno) {
  
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=mat, pathways=hallmark_pathways, phenotypes=as.factor(pheno)) 
  
  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)
  
  dd$metric <- "kendall_tau"
  dd$group1 <- "1"
  dd$group2 <- "2"
  return(dd)
}

splat_zeros <- function(mat, pheno) {
  
  #number of zeros
  g1 <- mat[,pheno == "1"]
  g2 <- mat[,pheno == "2"]
  
  g1_zero <- colSums(g1 == 0)
  g2_zero <- colSums(g2 == 0)
  
  g1box = data.frame(group = "1", value = g1_zero)
  g2box = data.frame(group = "2", value = g2_zero)
  
  zbox <- rbind(g1box, g2box)
}  
