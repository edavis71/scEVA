# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 3 - EVA breast tumor vs normal
# code developed by Emily Davis

###################################################
### code chunk: init
###################################################
sink = (file = "rout.log")
library(monocle)
library(GSReg)
load("hallmark_pathways.rda")
load("myeloid_pathways.rda")
load("pancancer_pathways.rda")
load("tcell_activation_pathways.rda")
load("breast_cancer_imputed_cdsBC1_8.rda")

###################################################
### code chunk 1: load biscuit imputed data BC 1-8
###################################################
tab <- read.csv('GSE114725_rna_imputed.csv',header=F, sep=",")

dat <- tab[1,]
genes <- t(dat[,-c(1,2,3,4,5)])
genes <- as.data.frame(genes)

# barcodes, cluster, patient for pdata
barcode <- tab[,c(5,1,2,3,4)]
barcode <- as.matrix(barcode) 
# assign first row to column names of matrix
colnames(barcode) <- barcode[1,]
# remove first row after setting it to col names of matrix
barcode <- barcode[-1,]
# assign first column to row names of matrix
row.names(barcode) <- barcode[,1]
barcode <- as.data.frame(barcode)

# remove first 5 rows
tab1 <- tab[,-c(1,2,3,4)]
tab1 <- as.matrix(tab1)
# assign first column to row names of matrix
row.names(tab1) <- tab1[,1]
# remove first column after setting it to row names of matrix
tab1 <- tab1[,-1]
# assign first row to column names of matrix
colnames(tab1) <- tab1[1,]
# remove first row after setting it to col names of matrix
tab1 <- tab1[-1,]
tab1 <- t(tab1)

# remove rownames for cds
rownames(tab1) <- c()
mat <- as.matrix(tab1)
matn <- apply(mat, 1, as.numeric) 
matn <- t(matn)
colnames(matn) <- colnames(mat)

###################################################
### code chunk 2: Create cell data set
###################################################

common_colnames <- c("gene_short_name")
fd <- new("AnnotatedDataFrame", data = genes)
colnames(fd) <- common_colnames
pd <- new("AnnotatedDataFrame", data = barcode)

cds = newCellDataSet(matn, featureData = fd, phenoData = pd)

###################################################
### code chunk 3: annotations
###################################################
# labeled from supp file
anno <- read.csv('bc1_8_anno.csv',header=T, sep=",")
common_colnames <- c("Cluster", "Celltype")
anno <- anno[c(1,6)]
colnames(anno) <- common_colnames

#check unique cell type clusters
unique(anno$Celltype)
anno$Cluster[(anno$Celltype == "T:CD4+NAIVE")]
anno$Cluster[(anno$Celltype == "T:CD8+NAIVE")]
anno$Cluster[(anno$Celltype == "NK:CD56+16+3+NKT")]
anno$Cluster[(anno$Celltype == "B:")]
anno$Cluster[(anno$Celltype == "NKT")]
anno$Cluster[(anno$Celltype == "T:CD8+EM")]
anno$Cluster[(anno$Celltype == "T:CD4+CM")]
anno$Cluster[(anno$Celltype == "mDC:")]
anno$Cluster[(anno$Celltype == "MACROPHAGE:")]
anno$Cluster[(anno$Celltype == "MONOCYTE:precursor")]
anno$Cluster[(anno$Celltype == "NEUTROPHIL:")]
anno$Cluster[(anno$Celltype == "NK:CD56+16+3-")]
anno$Cluster[(anno$Celltype == "MAST:")]
anno$Cluster[(anno$Celltype == "MONOCYTE:")]
anno$Cluster[(anno$Celltype == "pDC:")]
anno$Cluster[(anno$Celltype == "T:Reg")]
anno$Cluster[(anno$Celltype == "NK:CD56-16+3-")]
anno$Cluster[(anno$Celltype == "T:CD8+CM")]
anno$Cluster[(anno$Celltype == "T:CD4+EM")]


pData(cds)$Celltype == "NA"

pData(cds)$Celltype[(pData(cds)$cluster %in% c(1,2,17,19,20,22,60,81))] <- "T:CD4+NAIVE"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(3,5,43))] <- "T:CD8+NAIVE"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(4,12))] <- "NK:CD56+16+3+NKT"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(6,8,39,48,53,55,57,61,65,90))] <- "B:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(9))] <- "NKT"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(11,15,45,51,59,62,72,74,76,83,85))] <- "T:CD8+EM"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(16,47,66,75,95))] <- "T:CD4+CM"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(21,37,38))] <- "mDC:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(23,25,28))] <- "MACROPHAGE:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(30,35,63,78,88,92))] <- "NK:CD56+16+3-"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(31,44,73))] <- "MAST:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(36,40,54,67,68,84,86,91,94))] <- "MONOCYTE:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(41))] <- "pDC:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(46,56,80,82,87))] <- "T:Reg"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(41))] <- "pDC:"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(50))] <- "NK:CD56-16+3-"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(58,77))] <- "T:CD8+CM"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(70,71,93))] <- "T:CD4+EM"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(27))] <- "MONOCYTE:precursor"
pData(cds)$Celltype[(pData(cds)$cluster %in% c(28,89))] <- "NEUTROPHIL:"

# save cds with annotations
save(cds, file = "breast_cancer_imputed_cdsBC1_8.rda")

