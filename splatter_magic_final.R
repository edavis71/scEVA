# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 2 simulations
# code developed by Emily Davis

###################################################
### code chunk: init - run this script in terminal R
###################################################
library(splatter)
library(Rmagic)
library(ggplot2)
library(GSReg)
library(reshape2)
library(viridis)
library(GSBenchMark)

# load pathways
load("hallmark_pathways.rda")
# load list of genes
load("genelist.rda")
# load functions
source("splatter_functions.R")
set.seed(573763563)

###################################################
### code chunk: simulations: 2 identical groups, biased zeros
###################################################

# initate splatter parameters
params <- newSplatParams()
nGenes = 10000
params.dropout <- setParams(params, update = list(nGenes = nGenes, mean.rate = 0.2, mean.shape = 0.5, dropout.present = TRUE))
params <- setParams(params, update = list(nGenes = nGenes, mean.rate = 0.2, mean.shape = 0.5, dropout.present = FALSE))

# 1 group with and without dropout - no signal 
dropout <- splatSimulate(params.dropout)
nodropout <- splatSimulate(params)

dropout.matrix <- as.matrix(counts(dropout))
nod.matrix <- as.matrix(counts(nodropout))

combo <- cbind(dropout.matrix, nod.matrix)
phenotypes.combo <- rep(1:2, each = 100)
rownames(combo) <- genelist[1:10000]

# bias before imputing 
eva2 <- splat_eva(combo, phenotypes.combo)
b <- splat_zeros(combo, phenotypes.combo)

eva2$pvaluadj = p.adjust(eva2$pvalue, method='BH')
sig <- eva2[eva2$pvalue<0.05,]
common_cols <- c("estat1", "pathway", "group")
E1 <- sig[c(1,21,23)]
E2 <- sig[c(2,21,24)]
colnames(E1) <- common_cols
colnames(E2) <- common_cols
plotdat <- rbind(E1,E2)

pdf("final_spatter_identical_bias_100.pdf")
ggplot(plotdat, aes(x = reorder(pathway, estat1), y = estat1, color = group, shape=group, group=group)) +
  geom_point(alpha=0.6, size=3) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "pathway", y = "estat") + scale_color_manual(values=c("red3","dodgerblue4")) +
  geom_smooth(method = "lm", fill = NA)


pdf("splatter_dropout_100_zeros.pdf", width = 8, height = 5)
ggplot(b, aes(x=group, y=value, fill=group)) +  # This is the plot function
  geom_boxplot() + ggtitle("Number of Zeros Splatter") + 
  theme_classic() + scale_fill_manual(values=c("red3", "dodgerblue4")) +
  theme(legend.position="none") 
dev.off()

# impute with magic
premagic <- t(combo)
rownames(premagic) = make.names(rownames(premagic), unique=TRUE)
magicdat <- magic(premagic)

write.csv(magicdat, "final_splatter_biasedzeros_identical_100_magic.csv")
magicdat <- read.csv("final_splatter_biasedzeros_identical_100_magic.csv",header=T, sep=",")
rownames(magicdat) <- magicdat$X
magicdat <- magicdat[-c(1)]
postmagic <- t(magicdat) 

evam <- splat_eva(postmagic, phenotypes.combo)
c <- splat_zeros(postmagic, phenotypes.combo)

# no pathways significant after imputing
common_cols <- c("estat", "pathway", "group")
E1 <- evam[c(1,21,23)]
E2 <- evam[c(2,21,24)]
colnames(E1) <- common_cols
colnames(E2) <- common_cols
plotdat <- rbind(E1,E2)

pdf("final_spatter_identical_bias_100_magic.pdf")
ggplot(plotdat, aes(x = reorder(pathway, estat), y = estat, color = group, shape=group, group=group)) +
  geom_point(alpha=0.6, size=3) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "pathway", y = "estat") + scale_color_manual(values=c("red3","dodgerblue4")) +
  geom_smooth(method = "lm", fill = NA)
dev.off()

ggplot(c, aes(x=group, y=value, fill=group)) +  # This is the plot function
  geom_boxplot() + ggtitle("Number of Zeros Splatter") + theme_classic() + scale_fill_manual(values=c("red3", "dodgerblue4"))
dev.off()

###################################################
### code chunk: simulation: 2 identical groups, biased zeros, magic, random sampling
###################################################

# add randomized signal for pathways
for (i in 1:100) {
  for (j in 1:length(hallmark_pathways)) {
    include_list <- rownames(postmagic)[rownames(postmagic) %in% hallmark_pathways[[j]]]
    postmagic[include_list, i] <- postmagic[sample(include_list),i]
  }
}

eva_random <- splat_eva(postmagic, phenotypes.combo)
zeros_random <- splat_zeros(postmagic, phenotypes.combo)

eva_random$pvaluadj = p.adjust(eva_random$pvalue, method='BH')
sig <- eva_random[eva_random$pvaluadj<0.05,]
common_cols <- c("estat", "pathway", "group")
E1 <- sig[c(1,21,23)]
E2 <- sig[c(2,21,24)]
colnames(E1) <- common_cols
colnames(E2) <- common_cols
plotdat_random <- rbind(E1,E2)

dim(sig[sig$E1 > sig$E2,])

pdf("final_spatter_identical_bias_100_magic_randomized.pdf")
ggplot(plotdat_random, aes(x = reorder(pathway, estat),
                          y = estat, 
                          color = group, 
                          shape=group, 
                          group=group)) + 
  geom_point(alpha=0.6, size=3) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  labs(x = "pathway", y = "estat") + scale_color_manual(values=c("red3","dodgerblue4")) +
  geom_smooth(method = "lm", fill = NA)
dev.off()

pdf("splatter_dropout_100_randomized_magic_zeros.pdf", width = 8, height = 5)
ggplot(zeros_random, aes(x=group, y=value, fill=group)) +  # This is the plot function
  geom_boxplot() + ggtitle("Number of Zeros Splatter") + 
  theme_classic() + scale_fill_manual(values=c("red3", "dodgerblue4")) +
  theme(legend.position="none") 
dev.off()
