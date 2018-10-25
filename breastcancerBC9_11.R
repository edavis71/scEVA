# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 4 C EVA analysis between tumor t-cells
# code developed by Emily Davis

###################################################
### code chunk: init
###################################################
library(monocle)
library(GSReg)
library(reshape2)
library(ComplexHeatmap)
library(limma)
setwd("/Users/emilydavis/Desktop/breast_cancer/results")

# load data and pathways
load("breast_cancer_imputed_cdsBC9_11.rda")
load("myeloid_pathways.rda")

###################################################
### code chunk 1: load biscuit imputed data
###################################################
tab <- read.table(file = 'GSE114724_rna_imputed.tsv', sep = '\t', header = FALSE)

# barcodes, cluster, patient for pdata
barcode <- t(tab[c(2,1,3),])
# assign first row to column names of matrix
colnames(barcode) <- barcode[1,]
# remove first row after setting it to col names of matrix
barcode <- barcode[-1,]
# assign first column to row names of matrix
row.names(barcode) <- barcode[,1]
# remove first column after setting it to row names of matrix
barcode <- barcode[,-1]

# remove first and third row
tab1 <- tab[-c(1,3),]
# assign first column to row names of matrix
row.names(tab1) <- tab1[,1]
# remove first column after setting it to row names of matrix
tab1 <- tab1[,-1]
# assign first row to column names of matrix
colnames(tab1) <- rownames(barcode)
# remove first row after setting it to col names of matrix
tab1 <- tab1[-1,]

genes <- as.data.frame(tab[,1])
genes <- genes[-c(1,2,3),]
genes <- as.data.frame(genes)

barcode <- as.data.frame(barcode)

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
# labeled from distance heatmap fig 6
pData(cds)$Celltype == "NA"

pData(cds)$Celltype[(pData(cds)$Cluster %in% c(1,8,11,13,15,18,24,25,28,31,27))] <- "T:CD8+EM"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(10,19))] <- "T:CD8+CM"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(12,17,32,21))] <- "T:CD8+NAIVE"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(3,4))] <- "T:CD4+EM"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(5,6,9,26,30))] <- "T:CD4+CM"
pData(cds)$Celltype[which(pData(cds)$Cluster %in% c(7,16,22,33))] <- "T:CD4+NAIVE"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(2,14,20,29,34))] <- "T:Reg"
pData(cds)$Celltype[(pData(cds)$Cluster %in% c(23))] <- "NKT"

# save cds with annotations
save(cds, file = "breast_cancer_imputed_cdsBC9_11.rda")


###################################################
### code chunk 6: eva t cell by patient and celltype
###################################################

do_pathways2 <- function(pathway, sample1, sample2, celltype) {
  # subset cds by patient ID
  test<- cds[,(pData(cds)$Patient==sample1 | pData(cds)$Patient==sample2)]
  # subset cancer cells
  testc <- test[,(pData(test)$Celltype==celltype)]
  pData(testc)$compare <- ifelse(pData(testc)$Patient == sample1, 0,1)
  #celltype1 =0; celltype2 =1 in pData(testc)$compare
  # get expression matrix for pt cancer cells
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
  dd$patient1 <- sample1
  dd$patient2 <- sample2
  dd$celltype <- celltype
  return(dd)
}

a <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+EM")
b <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+EM")
c <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:CD4+EM")

d <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+NAIVE")
e <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+NAIVE")
f <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:CD4+NAIVE")

g <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+CM")
h <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+CM")
i <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:CD8+CM")

#stopped on k
j <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+EM")
k <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+EM")
l <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:CD8+EM")

m <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+NAIVE")
n <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+NAIVE")
o <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:CD8+NAIVE")

p <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "T:Reg")
q <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "T:Reg")
r <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "T:Reg")

s <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC9", sample2 = "BC10", celltype = "NKT")
t <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC10", sample2 = "BC11", celltype = "NKT")
u <- do_pathways2(pathway = hallmark_pathways, sample1 = "BC11", sample2 = "BC9", celltype = "NKT")


all <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)
all$pvaluadj <- p.adjust(all$pvalue,method='BH')
write.csv(all, file = "all_celltypes_eva_hallmark_BC9_11_results_atoi.csv")


all <- read.csv('all_celltypes_eva_myeloid_BC9_11_results.csv',header=T, sep=",")
all <- all[-c(1)]
#new <- all[!all$celltype == "NKT",]
common_cols <- c("estat", "pathway", "patient", "celltype")
all$col <- paste(all$patient1, all$celltype, sep="_")
E1 <- all[c(1,21,23,27)]
colnames(E1) <- common_cols

mat <- acast(E1, pathway ~ celltype , value.var='estat')

mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

# add heatmap annotation
df = data.frame(patient = rep(c("BC10", "BC11", "BC9"), each = 7),
                celltype = c("NKT", "T:CD4+EM", "T:CD4+NAIVE", "T:CD8+CM", "T:CD8+EM", "T:CD8+NAIVE",
                             "T:Reg", "NKT", "T:CD4+EM", "T:CD4+NAIVE", "T:CD8+CM", "T:CD8+EM", 
                             "T:CD8+NAIVE", "T:Reg", "NKT", "T:CD4+EM", "T:CD4+NAIVE", "T:CD8+CM",
                             "T:CD8+EM", "T:CD8+NAIVE", "T:Reg"))


# Green as highlight
greenfocus = c("#41AB5D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")

top = HeatmapAnnotation(df = df, col = list(patient = c("BC11" = "#619CFF", "BC9" = "#F8766D", "BC10" = "#00BA38"), 
                                            celltype = c("NKT" = "#252525", 
                                                         "T:CD4+EM" = "#525252", 
                                                         "T:CD4+NAIVE" = "#737373",
                                                         "T:CD8+CM" = "#969696", 
                                                         "T:CD8+EM" = "#BDBDBD", 
                                                         "T:CD8+NAIVE" = "#D9D9D9", 
                                                         "T:Reg" = "#F0F0F0")))

pdf("celltypes_eva_myeloid_BC9_11_legend_complexheatmap_scale.pdf",width=15,height=10, paper = "USr")
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"), top_annotation = top, row_names_side = "left", 
        row_dend_side = "right", width = unit(80, "mm"), row_names_gp = gpar(fontsize = 8))
        #column_title = "T cell dysregultation by pathway: myeloid")
dev.off()


###################################################
### code chunk 7: DE t cell by patient and celltype
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
  
  #testc <- estimateSizeFactors(testc)
  #testc <- estimateDispersions(testc)
  
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

gsea <- function(pathways, data){
   
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
  return(stats)
}

barcode_plot <- function(stats){
  
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

a <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+EM")
a_gst <- gsea(pathways = myeloid, data = a)
barcode_plot(a_gst)

b <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+EM")
b_gst <- gsea(pathways = myeloid, data = b)
barcode_plot(b_gst)

c <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD4+EM")
c_gst <- gsea(pathways = myeloid, data = c)
barcode_plot(c_gst)

d <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD4+NAIVE")
d_gst <- gsea(pathways = myeloid, data = d)
barcode_plot(d_gst)

e <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD4+NAIVE")
e_gst <- gsea(pathways = myeloid, data = e)
barcode_plot(e_gst)

f <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD4+NAIVE")
f_gst <- gsea(pathways = myeloid, data = f)
barcode_plot(f_gst)

g <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+CM")
g_gst <- gsea(pathways = myeloid, data = g)
barcode_plot(g_gst)

h <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+CM")
h_gst <- gsea(pathways = myeloid, data = h)
barcode_plot(h_gst)

i <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+CM")
i_gst <- gsea(pathways = myeloid, data = i)
barcode_plot(i_gst)

j <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+EM")
j_gst <- gsea(pathways = myeloid, data = j)
barcode_plot(j_gst)

k <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+EM")
k_gst <- gsea(pathways = myeloid, data = k)
barcode_plot(k_gst)

l <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+EM")
l_gst <- gsea(pathways = myeloid, data = l)
barcode_plot(l_gst)

m <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:CD8+NAIVE")
m_gst <- gsea(pathways = myeloid, data = m)
barcode_plot(m_gst)

n <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:CD8+NAIVE")
n_gst <- gsea(pathways = myeloid, data = n)
barcode_plot(n_gst)

o <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:CD8+NAIVE")
o_gst <- gsea(pathways = myeloid, data = o)
barcode_plot(o_gst)

p <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "T:Reg")
p_gst <- gsea(pathways = myeloid, data = p)
barcode_plot(p_gst)

q <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "T:Reg")
q_gst <- gsea(pathways = myeloid, data = q)
barcode_plot(q_gst)

r <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "T:Reg")
r_gst <- gsea(pathways = myeloid, data = r)
barcode_plot(r_gst)

s <- de_test(sample1 = "BC9", sample2 = "BC10", celltype = "NKT")
s_gst <- gsea(pathways = myeloid, data = s)
barcode_plot(s_gst)

t <- de_test(sample1 = "BC10", sample2 = "BC11", celltype = "NKT")
t_gst <- gsea(pathways = myeloid, data = t)
barcode_plot(t_gst)

u <- de_test(sample1 = "BC9", sample2 = "BC11", celltype = "NKT")
u_gst <- gsea(pathways = myeloid, data = u)
barcode_plot(u_gst)

# "T:CD4+EM", "T:CD4+NAIVE", "T:CD8+CM", "T:CD8+EM", "T:CD8+NAIVE", "T:Reg", "NKT"

# write gsea data
gsea_dat_down <- rbind(a_gst, b_gst, c_gst, d_gst, e_gst, f_gst, g_gst, h_gst, i_gst, j_gst, k_gst, l_gst, 
      m_gst, n_gst, o_gst, p_gst, q_gst, r_gst, s_gst, t_gst, u_gst)

#write.csv(gsea_dat, file ="all_BC9_11_GSEA_myeloid.csv")


# either
sig_gsea <- gsea_dat[gsea_dat$gst < 0.05,]
# up
sig_gsea_up <- gsea_dat_up[gsea_dat_up$gst < 0.05,]
# down
sig_gsea_down <- gsea_dat_down[gsea_dat_down$gst < 0.05,]

library(gridExtra)
pdf("all_BC9_11_sig_GSEA_myeloid_results", width=10,height=38)
grid.table(sig_gsea)
dev.off()

# compare lists
library(gplots)

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
library(dplyr)
dim(inner_join(eva, sig_gsea))
dim(inner_join(sig_gsea_up, sig_gsea_down))

eva_col <- paste(eva$pathway_name, eva$patient1, eva$patient2, eva$celltype, sep="_")
gsea_col_up <- paste(sig_gsea_up$pathway_name, sig_gsea_up$patient1, sig_gsea_up$patient2, sig_gsea_up$celltype, sep="_")
gsea_col_down <- paste(sig_gsea_down$pathway_name, sig_gsea_down$patient1, sig_gsea_down$patient2, sig_gsea_down$celltype, sep="_")

library(Vennerable)
v <- list("GSEA UP"=gsea_col_up,"GSEA DOWN"=gsea_col_down,"EVA"=eva_col)
VennRaw <- Venn(v)

pdf("venn_BC9_11_EVA_GSEA_myeloid_results_up_down", width=5,height=5)
plot(VennRaw, doWeights = FALSE)
dev.off()
