# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 1 overview matrix, tsne, barplot
# code developed by Emily Davis

###################################################
### code chunk 1: r init
###################################################
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dplyr)
library(magrittr)
library(ggrepel)
library(Rtsne)

###################################################
### code chunk: overview fig matrix
###################################################

N = 50
M = 50
matvar <- matrix( rnorm(N*M,mean=0,sd=1), N, M) 

# remove first column - adjust sd to change variation level
ma <- matvar[,1]
reps <- t(do.call(rbind, replicate(25, ma, simplify=FALSE)))

# variable data
jitv <- reps + matrix(rnorm(nrow(reps) * ncol(reps), sd = 50), ncol = ncol(reps))
scalevar <- apply(jitv,2,scale)
meltvar <- melt(scalevar)

var <- ggplot(meltvar, aes(Var2,Var1)) + geom_tile(aes(fill = (value)), colour = "white") + 
  scale_fill_gradient2(low = "steelblue3", mid = "khaki1", high = "red1") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  ggtitle("variable") +labs(x = "cells", y = "pathway genes")

# data with little variation 
jitc <- reps + matrix(rnorm(nrow(reps) * ncol(reps), sd = 0.25), ncol = ncol(reps))
scalec <- apply(jitc,2,scale)
meltc <- melt(scalec)

const <- ggplot(meltc, aes(Var2,Var1)) +
  geom_tile(aes(fill = (value)), colour = "white") + 
  scale_fill_gradient2(low = "steelblue3", mid = "khaki1", high = "red1") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  ggtitle("non-variable") +labs(x = "cells", y = "pathway genes")

# plot both matrices
grid.arrange(const, var, ncol = 2)

# combine into one matrix
mat <- cbind(jitc, jitv)
dim(mat)
colnames(mat) <- rep(1:50, each = 1)
mat <- mat[,-c(1:5)]
scaled <- apply(mat,2,scale)
meltd <- melt(scaled)

pdf("overview_fig_all_matrix.pdf",width=15,height=10, paper = "USr")
ggplot(meltd, aes(Var2,Var1)) +
  geom_tile(aes(fill = (value)), colour = "white") + 
  scale_fill_gradient2(low = "steelblue3", mid = "khaki1", high = "red1") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(x = "samples", y = "pathway genes")

dev.off()

# generate top legend
df = as.data.frame(rep(c("a", "b"), each = 25)) 
colnames(df) <- c("group")
ha = HeatmapAnnotation(df = df, col = list(group = c("a" =  "blue", "b" = "red")))
draw(ha, 1:50)

###################################################
### code chunk: overview tsne 
###################################################

set.seed(123456)
N = 2000
D = 200
data.norm = matrix(rnorm(N * D, 2), N)
groups.probs = runif(2)
groups = sample(1:2, N, TRUE, groups.probs/sum(groups.probs))
for (gp in unique(groups)) {
  dev = rep(1, D)
  dev[sample.int(D, 3)] = runif(3, -10, 10)
  data.norm[which(groups == gp), ] = data.norm[which(groups == 
                                                       gp), ] %*% diag(dev)
}
info.norm = tibble(truth = factor(groups))

pca.norm = prcomp(data.norm)
info.norm %<>% cbind(pca.norm$x[, 1:4])

tsne.norm = Rtsne(pca.norm$x, pca = FALSE)
info.norm %<>% mutate(tsne1 = tsne.norm$Y[, 1], tsne2 = tsne.norm$Y[,2])

pdf("overview_fig_tsne.pdf",width=6.5,height=5)
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = truth)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  scale_color_manual(breaks = c("1", "2"),
                     values=c("blue", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

###################################################
### code chunk: overview bar plot of estats
###################################################

phenotype <- c("1", "2")
estat <- c(0.10, 0.22)
df <- data.frame(phenotype, estat)

pdf("overview_fig_bar.pdf",width=5,height=5)
ggplot(data = df, aes(x=phenotype, y=estat, fill=phenotype)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values =c("blue", "red")) +
  geom_signif(y_position=c(0.25), xmin=c(1), xmax=c(2),annotation=c("***"), tip_length =0.2)
dev.off()
