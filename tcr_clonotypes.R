# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 4 A and B, tcr clonality, richness, and morisita matrix
# code developed by Emily Davis

###################################################
### code chunk: init
###################################################

library(tcrSeqR)
library(reshape2)
library(ggplot2)
library(viridis)
setwd("/Users/emilydavis/Desktop/breast_cancer/GSE114724_RAW")

###################################################
### code chunk 1: prep files
###################################################

file1 <- read.csv('GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv',header=T, sep=",")
file1$sample <- "BC9_01"
file2 <- read.csv('GSM3148581_BC09_TUMOR2_filtered_contig_annotations.csv',header=T, sep=",")
file2$sample <- "BC9_02"
file3 <- read.csv('GSM3148582_BC10_TUMOR1_filtered_contig_annotations.csv',header=T, sep=",")
file3$sample <- "BC10_01"
file4 <- read.csv('GSM3148583_BC11_TUMOR1_filtered_contig_annotations.csv',header=T, sep=",")
file4$sample <- "BC11_01"
file5 <- read.csv('GSM3148584_BC11_TUMOR2_filtered_contig_annotations.csv',header=T, sep=",")
file5$sample <- "BC11_02"

# get count data in iseq format
f1 <- file1[c(14,13,15,19)]
f2 <- file2[c(14,13,15,19)]
f3 <- file3[c(14,13,15,19)]
f4 <- file4[c(14,13,15,19)]
f5 <- file5[c(14,13,15,19)]

all <- rbind(f1,f2,f3,f4,f5)
common_colnames <- c("nt", "aa", "reads", "sample")
colnames(all) <- common_colnames

aggdata <-aggregate(all$reads, by=list(all$nt,all$aa,all$sample), 
                    FUN=mean, na.rm=TRUE)
common_colnames <- c("nt", "aa", "sample", "reads")
colnames(aggdata) <- common_colnames

BC9_01 <- aggdata[aggdata$sample == "BC9_01",]
BC9_02 <- aggdata[aggdata$sample == "BC9_02",]
BC10_01 <- aggdata[aggdata$sample == "BC10_01",]
BC11_01 <- aggdata[aggdata$sample == "BC11_01",]
BC11_02 <- aggdata[aggdata$fsample == "BC11_02",]

all_nt <- aggdata$nt
all_aa <- aggdata$aa

list <- data.frame(nt=all_nt,aa=all_aa)
dim(list)
list <- unique(list)
dim(list)

dat1 <- merge(list, BC9_01[, c("nt", "reads")], by="nt", all.x=TRUE)
dat1$reads[is.na(dat1$reads)] <- 0
dat2 <- merge(list, BC9_02[, c("nt", "reads")], by="nt", all.x=TRUE)
dat2$reads[is.na(dat2$reads)] <- 0
dat3 <- merge(list, BC10_01[, c("nt", "reads")], by="nt", all.x=TRUE)
dat3$reads[is.na(dat3$reads)] <- 0
dat4 <- merge(list, BC11_01[, c("nt", "reads")], by="nt", all.x=TRUE)
dat4$reads[is.na(dat4$reads)] <- 0
dat5 <- merge(list, BC11_02[, c("nt", "reads")], by="nt", all.x=TRUE)
dat5$reads[is.na(dat5$reads)] <- 0

list$BC9_01 <- dat1[,(3)]
list$BC9_02 <- dat2[,(3)]
list$BC10_01 <- dat3[,(3)]
list$BC11_01 <- dat4[,(3)]
list$BC11_02 <- dat5[,(3)]

# remove none line
list <- list[!list$nt == "None",]

#write.csv(list, file = "summarized_tcell_reads.csv")

# metadata dictionary
fn <- c("BC9_01", "BC9_02", "BC10_01", "BC11_01", "BC11_02")
patient <- c("BC9", "BC9", "BC10", "BC11", "BC11")
type <- patient
meta <- data.frame(fn = fn, patient = patient, type=type)
meta[] <- lapply(meta, as.character)

###################################################
### code chunk 2: tcrseqr clonality and richness
###################################################

# remove reads with stop codons 
ds_nt <- iseqr_make_tcr(list,meta, remove=T)

# aggregate the data - this collapses synonymous nucleotide sequences
ds_agg <- iseqr_aggregate(ds_nt)

# Clonality for one sample
clonality(assay(ds_agg)[,1])

#Clonality for all samples
clonality(ds_agg, merge=F)

# and these can be merged directly into the object's metadata
ds_agg <- clonality(ds_agg)
ds_agg <- richness(ds_agg)
ds_agg <- total(ds_agg)
ds_agg$Clonality
ds_agg$Richness

clone <- data.frame(patient = ds_agg$patient, clonality=ds_agg$Clonality, richness = ds_agg$Richness)

clone$patient <- factor(clone$patient, levels=c("BC9", "BC10", "BC11"))

plot1 <- ggplot(clone, aes(x = reorder(patient, clonality), y = clonality, color = patient)) +
  geom_point(alpha=0.6,size=3) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  #stat_summary(fun.data = "mean_cl_boot", size = 1.0, shape = 3) +
  ggtitle("Tumor T Cell Clonality") + labs(x = "patient")


plot2 <- ggplot(clone, aes(x = reorder(patient, clonality), y = richness, color = patient)) +
  geom_point(alpha=0.6, size=3) + 
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  #stat_summary(fun.data = "mean_cl_boot", size = 1.0, shape = 3) +
  ggtitle("Tumor T Cell Richness") + labs(x = "patient")

require(gridExtra)
pdf("BC9_11_clonality_richness_ordered_grid.pdf",width=8,height=5)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

###################################################
### code chunk 3: morisita similarity matrix
###################################################

# morisita matrix of replicates
morisita(assay(ds_agg)[,1], assay(ds_agg)[,3])

# build empty matrix
mat <- matrix(nrow = 5, ncol = 5)
colnames(mat) = colnames(assay(ds_agg))
rownames(mat) = colnames(assay(ds_agg))

  for(k in 1:ncol(assay(ds_agg)))
  {
    for (l in 1:ncol(assay(ds_agg)))
    {
     # dist = cenken(V[,k], L[,k], V[,l], L[,l])
      dist = morisita(assay(ds_agg)[,k], assay(ds_agg)[,l])
  
      #colname of V[k]
      sample1 <- colnames(assay(ds_agg))[[k]][1]
      #colname of V[l]
      sample2 <- colnames(assay(ds_agg))[[l]][1]
      
      mat[sample1, sample2]  = dist
    }
  }    

# make data in long format for ggplot
melted_mat <- melt(mat)
# order samples numerically
melted_mat$X1 <- factor(melted_mat$X1 , levels=c("BC9_01","BC9_02","BC10_01","BC11_01","BC11_02"))
melted_mat$X2 <- factor(melted_mat$X2 , levels=c("BC9_01","BC9_02","BC10_01","BC11_01","BC11_02"))

ggplot(data = melted_mat, aes(x=X1, y=X2, fill=log2(value))) + 
  geom_tile() + scale_fill_viridis(option="plasma") 
