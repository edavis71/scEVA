# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 4 C, GSEA analysis
# code developed by Emily Davis

###################################################
### code chunk: init
###################################################

library(monocle)
library(reshape2)
library(limma)
library(gridExtra)
library(Vennerable)
library(dplyr)
library(gplots)
setwd("/Users/emilydavis/Desktop/breast_cancer")
setwd("/Users/emilydavis/Desktop/breast_cancer/results")

###################################################
### code chunk 2: Create cell data set
###################################################

load("breast_cancer_imputed_cdsBC9_11.rda")

# need cds in neg binomial for DE
exprsdat <- exp(exprs(cds))
rownames(exprsdat) <- fData(cds)$gene_short_name
fdat <- fData(cds)
rownames(fdat) <- fData(cds)$gene_short_name
pdat <- pData(cds)

common_colnames <- c("gene_short_name", "num_cells_expressed")
fd <- new("AnnotatedDataFrame", data = fdat)
colnames(fd) <- common_colnames
pd <- new("AnnotatedDataFrame", data = pdat)

cds = newCellDataSet(matn, featureData = fd, phenoData = pd)

cds <- newCellDataSet(exprsdat,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())

# save neg binomial cds
save(cds, file = "breast_cancer_imputed_cdsBC9_11_gsea.rda")

###################################################
### code chunk: DE functions
###################################################

de_test <- function(sample1, sample2, celltype) {
  
  # subset for 2 patients and single celltype
  test<- cds[,(pData(cds)$Patient==sample1 | pData(cds)$Patient==sample2)]
  testc <- test[,(pData(test)$Celltype==celltype)]
  
  #Quantifying gene expression by patient, annotating expression changes; lower pt # first
  testc_sample1_count <- sum(pData(testc)$Patient==sample1 , na.rm=TRUE)
  testc_sample2_count <- sum(pData(testc)$Patient==sample2, na.rm=TRUE)
  fData(testc)$Total_mRNAs_per_gene_testc <- Matrix::rowSums(exprs(testc))
  fData(testc)$Average_mRNAs_testc_sample1 <- Matrix::rowSums(exprs(testc[,pData(testc)$Patient==sample1]))/testc_sample1_count
  fData(testc)$Average_mRNAs_testc_sample2 <- Matrix::rowSums(exprs(testc[,pData(testc)$Patient==sample2]))/testc_sample2_count
  fData(testc)$Expression_change_testc <- ifelse(fData(testc)$Average_mRNAs_testc_sample1 < fData(testc)$Average_mRNAs_testc_sample2, "Upregulated", ifelse(fData(testc)$Average_mRNAs_testc_sample1 > fData(testc)$Average_mRNAs_testc_sample2, "Downregulated", "Unchanged"))
  fData(testc)$Fold_change_testc <- fData(testc)$Average_mRNAs_testc_sample2/fData(testc)$Average_mRNAs_testc_sample1
  
  testc <- estimateSizeFactors(testc)
  testc <- estimateDispersions(testc)
  
  testc<-detectGenes(testc,min_expr=1)
  
  numCellThreshold<-5 
  testc.expressed_genes<-row.names(subset(fData(testc),num_cells_expressed >= numCellThreshold))
  head(testc.expressed_genes)
  
  DE_test_res <- differentialGeneTest(testc[testc.expressed_genes,], 
                                      fullModelFormulaStr = "~num_genes_expressed+Patient", 
                                      reducedModelFormulaStr = "~num_genes_expressed",
                                      cores=detectCores())
  
  DE_test_sigGenes <- subset(DE_test_res, qval < 0.05)
  
  DE_test_sigGenes$sample1 = sample1
  DE_test_sigGenes$sample2 = sample2
  DE_test_sigGenes$celltype = celltype
  return(DE_test_sigGenes)
}

gsea_up <- function(pathways, data){
  
  stats = NULL
  for (i in 1:length(pathways)) {
    pathway = pathways[i]
    pathway_name = names(pathway)
    pathway_genes = pathway[[1]]  
    
    path <- ifelse(data$gene_short_name %in% pathway_genes, TRUE,FALSE)
    data$logfc <- log2(data$Fold_change_testc)
    logfc <- data[c(15)]
    logfc <- as.matrix(logfc)
    logfc <- as.numeric(logfc)
    
    gst <- wilcoxGST(index=path, statistics=logfc, alternative='up')
    stats = rbind(stats, data.frame(gst, pathway_name))
  }
  
  stats$sample1 <- data$sample1[1]
  stats$sample2 <- data$sample2[1]
  stats$celltype <- data$celltype[1]
  stats$gsea <- "up"
  return(stats)
}

gsea_down <- function(pathways, data){
  
  stats = NULL
  for (i in 1:length(pathways)) {
    pathway = pathways[i]
    pathway_name = names(pathway)
    pathway_genes = pathway[[1]]  
    
    path <- ifelse(data$gene_short_name %in% pathway_genes, TRUE,FALSE)
    data$logfc <- log2(data$Fold_change_testc)
    logfc <- data[c(15)]
    logfc <- as.matrix(logfc)
    logfc <- as.numeric(logfc)
    
    gst <- wilcoxGST(index=path, statistics=logfc, alternative='down')
    stats = rbind(stats, data.frame(gst, pathway_name))
  }
  
  stats$sample1 <- data$sample1[1]
  stats$sample2 <- data$sample2[1]
  stats$celltype <- data$celltype[1]
  stats$gsea <- "down"
  return(stats)
}

barcode_plot <- function(stats, data){
  
  # significant gsea pathways
  sig <- stats[stats$gst < 0.05,]
  
  pdfname <- paste("barcode_plots", sig$sample1, "vs", sig$sample2, sig$celltype, sep="_")
  pdf(pdfname, width=15,height=8, paper = "USr")
  for (i in 1:length(sig)) {
    pathway_genes = myeloid[[sig$pathway_name[i]]]
    pathway_name = sig$pathway_name[i]
    
    path <- ifelse(data$gene_short_name %in% pathway_genes, TRUE,FALSE)
    data$logfc <- log2(data$Fold_change_testc)
    logfc <- data[c(15)]
    logfc <- as.matrix(logfc)
    logfc <- as.numeric(logfc)
    
    plotname <- paste("logFC:", pathway_name, sig$sample1[1], "vs", sig$sample2[1], sig$celltype[1], sep="_")
    
    barcodeplot(index=path, statistics=logfc, main = plotname)
  }
  dev.off()
}

###################################################
### code chunk: DE tests
###################################################
load("breast_cancer_imputed_cdsBC9_11_gsea.rda")

a <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+EM")
a_up <- gsea_up(pathways = myeloid, data = a)
a_down <- gsea_down(pathways = myeloid, data = a)
barcode_plot(stats = a_down, data = a)

b <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+EM")
b_up <- gsea_up(pathways = myeloid, data = b)
b_down <- gsea_down(pathways = myeloid, data = b)

c <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD4+EM")
c_up <- gsea_up(pathways = myeloid, data = c)
c_down <- gsea_down(pathways = myeloid, data = c)
barcode_plot(stats = c_up, data = c)


d <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+NAIVE")
d_up <- gsea_up(pathways = myeloid, data = d)
d_down <- gsea_down(pathways = myeloid, data = d)

e <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+NAIVE")
e_up <- gsea_up(pathways = myeloid, data = e)
e_down <- gsea_down(pathways = myeloid, data = e)

f <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD4+NAIVE")
f_up <- gsea_up(pathways = myeloid, data = f)
f_down <- gsea_down(pathways = myeloid, data = f)

g <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+CM")
g_up <- gsea_up(pathways = myeloid, data = g)
g_down <- gsea_down(pathways = myeloid, data = g)

h <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+CM")
h_up <- gsea_up(pathways = myeloid, data = h)
h_down <- gsea_down(pathways = myeloid, data = h)

i <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+CM")
i_up <- gsea_up(pathways = myeloid, data = i)
i_down <- gsea_down(pathways = myeloid, data = i)

j <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+EM")
j_up <- gsea_up(pathways = myeloid, data = j)
j_down <- gsea_down(pathways = myeloid, data = j)

k <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+EM")
k_up <- gsea_up(pathways = myeloid, data = k)
k_down <- gsea_down(pathways = myeloid, data = k)

l <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+EM")
l_up <- gsea_up(pathways = myeloid, data = l)
l_down <- gsea_down(pathways = myeloid, data = l)

m <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+NAIVE")
m_up <- gsea_up(pathways = myeloid, data = m)
m_down <- gsea_down(pathways = myeloid, data = m)

n <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+NAIVE")
n_up <- gsea_up(pathways = myeloid, data = n)
n_down <- gsea_down(pathways = myeloid, data = n)

o <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+NAIVE")
o_up <- gsea_up(pathways = myeloid, data = o)
o_down <- gsea_down(pathways = myeloid, data = o)

p <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:Reg")
p_up <- gsea_up(pathways = myeloid, data = p)
p_down <- gsea_down(pathways = myeloid, data = p)

q <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:Reg")
q_up <- gsea_up(pathways = myeloid, data = q)
q_down <- gsea_down(pathways = myeloid, data = q)

r <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:Reg")
r_up <- gsea_up(pathways = myeloid, data = r)
r_down <- gsea_down(pathways = myeloid, data = r)

s <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "NKT")
s_up <- gsea_up(pathways = myeloid, data = s)
s_down <- gsea_down(pathways = myeloid, data = s)

t <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "NKT")
t_up <- gsea_up(pathways = myeloid, data = t)
t_down <- gsea_down(pathways = myeloid, data = t)

u <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "NKT")
u_up <- gsea_up(pathways = myeloid, data = u)
u_down <- gsea_down(pathways = myeloid, data = u)

# write gsea data
gsea_dat_up <- rbind(a_up, b_up, c_up, d_up, e_up, f_up, g_up, h_up, i_up, j_up, k_up, l_up, 
                       m_up, n_up, o_up, p_up, q_up, r_up, s_up, t_up, u_up)

gsea_dat_down <- rbind(a_down, b_down, c_down, d_down, e_down, f_down, g_down, h_down, i_down, j_down, k_down, l_down, 
                       m_down, n_down, o_down, p_down, q_down, r_down, s_down, t_down, u_down)

write.csv(gsea_dat_up, file ="all_BC9_11_GSEA_myeloid_up.csv")
write.csv(gsea_dat_down, file ="all_BC9_11_GSEA_myeloid_down.csv")

# up
sig_gsea_up <- gsea_dat_up[gsea_dat_up$gst < 0.05,]
# down
sig_gsea_down <- gsea_dat_down[gsea_dat_down$gst < 0.05,]

# read in eva data
all <- read.csv('all_celltypes_eva_myeloid_BC9_11_results.csv',header=T, sep=",")
all <- all[-c(1)]

eva_sig <- all[all$pvaluadj < 0.05,]
switch <- eva_sig[which(eva_sig$patient1 == "BC11" & eva_sig$patient2 == "BC9"),]
switch$patient1 = "BC9"
switch$patient2 = "BC11"

sig_eva <- eva_sig[!(eva_sig$patient1 == "BC11" & eva_sig$patient2 == "BC9"),]
eva <- rbind(sig_eva, switch)


eva$col <- paste(eva$pathway_name, eva$patient1, eva$patient2, eva$celltype, sep="_")
sig_gsea_up$col <- paste(sig_gsea_up$pathway_name, sig_gsea_up$sample1, sig_gsea_up$sample2, sig_gsea_up$celltype, sep="_")
sig_gsea_down$col <- paste(sig_gsea_down$pathway_name, sig_gsea_down$sample1, sig_gsea_down$sample2, sig_gsea_down$celltype, sep="_")

v <- list("GSEA UP"=gsea_col_up,"GSEA DOWN"=gsea_col_down,"EVA"=eva_col)
VennRaw <- Venn(v)

# all pathways and celltypes
pdf("venn_BC9_11_EVA_GSEA_myeloid_results_all", width=5,height=5)
plot(VennRaw, doWeights = FALSE)
dev.off()

# table for supplement
pdf("all_BC9_11_sig_down_GSEA_myeloid_results", width=10,height=38)
grid.table(sig_gsea_down)
dev.off()

# separate comparisons and cell types

sub_e <- eva[eva$patient1 == "BC9" & eva$patient2 == "BC11",]
sub_up <- sig_gsea_up[sig_gsea_up$sample1 == "BC9" & sig_gsea_up$sample2 == "BC11",]
sub_down <- sig_gsea_down[sig_gsea_down$sample1 == "BC9" & sig_gsea_down$sample2 == "BC11",]

sub_e <- sub_e[sub_e$celltype == "T:CD8+EM",]
sub_up <- sub_up[sub_up$celltype == "T:CD8+EM",]
sub_down <- sub_down[sub_down$celltype == "T:CD8+EM",]

v <- list("GSEA UP"=sub_up$col,"GSEA DOWN"=sub_down$col,"EVA"=sub_e$col)
VennRaw <- Venn(v)

pdf("venn_BC9_11_EVA_GSEA_myeloid_results_all", width=5,height=5)
plot(VennRaw, doWeights = FALSE)
dev.off()


# compare lists
sig_eva <- all[all$pvaluadj < 0.05,]
sig_eva <- sig_eva[c(22,24,25,26)]

sig_gsea_up <- sig_gsea_up[-c(1)]
sig_gsea_down <- sig_gsea_down[-c(1)]
common_colnames <- c("pathway_name", "patient1", "patient2", "celltype")
colnames(sig_gsea_up) <- common_colnames
colnames(sig_gsea_down) <- common_colnames

switch <- sig_eva[which(sig_eva$patient1 == "BC11" & sig_eva$patient2 == "BC9"),]
switch$patient1 = "BC9"
switch$patient2 = "BC11"

sig_eva <- sig_eva[!(sig_eva$patient1 == "BC11" & sig_eva$patient2 == "BC9"),]
eva <- rbind(sig_eva, switch)

# make venn diagram 
dim(inner_join(eva, sig_gsea))
dim(inner_join(sig_gsea_up, sig_gsea_down))

eva_col <- paste(eva$pathway_name, eva$patient1, eva$patient2, eva$celltype, sep="_")
gsea_col_up <- paste(sig_gsea_up$pathway_name, sig_gsea_up$patient1, sig_gsea_up$patient2, sig_gsea_up$celltype, sep="_")
gsea_col_down <- paste(sig_gsea_down$pathway_name, sig_gsea_down$patient1, sig_gsea_down$patient2, sig_gsea_down$celltype, sep="_")


v <- list("GSEA UP"=gsea_col_up,"GSEA DOWN"=gsea_col_down,"EVA"=eva_col)
VennRaw <- Venn(v)

pdf("venn_BC9_11_EVA_GSEA_myeloid_results_up_down", width=5,height=5)
plot(VennRaw, doWeights = FALSE)
dev.off()
