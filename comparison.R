
###################################################
### code chunk: init
###################################################
library(splatter)
library(Rmagic)
library(scde)
library(GSReg)

###################################################
### code chunk: functions
###################################################

splat_eva <- function(mat, pheno, path) {
  
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=mat, pathways=path, phenotypes=as.factor(pheno)) 
  
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

###################################################
### code chunk: simulated data
###################################################

source('differentialVarianceSimulation_R.R')
# generate pathways
genes_A <- paste0('gene.A.', seq_len(30))
genes_B <- paste0('gene.B.', seq_len(40))
genes_C <- paste0('gene.C.', seq_len(35))
nCells <- 500
# set pathway variance
pwyA <- Pathway(genes_A, treatment_variance=1)
pwyB <- Pathway(genes_B, treatment_variance=1)
pwyC <- Pathway(genes_C, treatment_variance=1)
allPathways <- list(pwyA, pwyB, pwyC)
# simulate data with splatter
sim <- differentialVarianceSimulation(allPathways, nCellsPerGroup=nCells)

###################################################
### code chunk: imputation and EVA
###################################################

# get counts from splatter object 
counts <- BiocGenerics::counts(sim) 
phenotypes <- rep(1:2, each = 500)
# imputes scRNA count data
magic_counts <- magic(t(counts))
magic_counts <- t(magic_counts[["result"]])

# generate list of pathways 
myPathways <- list(genes_A, genes_B, genes_C)
names(myPathways) <- c("pwyA", "pwyB", "pwyC")

# EVA for differential variation
# variation present in pwyA
eva <- splat_eva(magic_counts, phenotypes, path = myPathways)
eva$pvalue < 0.05

###################################################
### code chunk: CV GSEA 
###################################################

# calculate coefficient of variation by gene
# then perfrom GSEA 
group1 <- counts
sd1 <- rowSds(group1)
mean1 <- rowMeans(group1)
cv1 <- sd1/mean1

gsea_cv <- data.frame()
for (i in 1:length(myPathways)){
  p <- geneSetTest(index = which(names(cv1) %in% myPathways[[i]]), statistics = cv1, alternative = "mixed")
  path <- names(myPathways)[i]
  gsea_cv[i, 1] <- path
  gsea_cv[i, 2] <- p
}
gsea_cv

###################################################
### code chunk: PAGODA
###################################################

# PAGODA analysis for pathway overdispersion
# set min.size.entries to number of total genes 
knn <- knn.error.models(counts, k = ncol(counts)/2, n.cores = 1, min.count.threshold = 1, min.nonfailed = 1, min.size.entries = 105, max.model.plots = 10)
varinfo <- pagoda.varnorm(knn, counts = counts, trim = 1/ncol(counts), max.adj.var = 1, n.cores = 1, plot = TRUE)
# convert pathway list to environment
path.env <- list2env(myPathways)
pwpca <- pagoda.pathway.wPCA(varinfo, path.env, n.components = 1, n.cores = 1)
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
2*pnorm(-abs(df$z))

