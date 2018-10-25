# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# supplemental figures
# code developed by Emily Davis

###################################################
### code chunk 1: r init
###################################################
library(ggplot2)
library(reshape2)
library(GSBenchMark)
library(GSReg)
library(philentropy)
setwd("/Users/emilydavis/Desktop/alldist")

################################################# 
### code chunk: load data
################################################# 
data(GSBenchMarkDatasets) 
DataSetStudy=GSBenchMark.Dataset.names[[9]]  #select the 9th dataset, squamous_GDS2520
data(list=DataSetStudy)
genenames=rownames(exprsdata)
print(phenotypesLevels)

load('hallmark_pathways.rda')

################################################# 
### code chunk: functions
################################################# 

# list of distance metrics to use 
dlist <- getDistMethods()
remove = c("minkowski", "non-intersection", "kulczynski_s", 
           "bhattacharyya", "hellinger", "matusita", "kullback-leibler", 
           "k_divergence", "motyka", "additive_symm", "topsoe", 
           "jensen_difference")
dlist <- dlist[!dlist %in% remove]
dlist <- append(dlist, "kendall_tau")

# function to calculate eva statistics
eva_stats <- function(df, dist, pheno) {
  phenotypes = as.factor(pheno)
  samplesC1 <- colnames(dist[,phenotypes==0])
  samplesC2 <- colnames(dist[,phenotypes==1])
  
  dist1 <- dist[samplesC1,samplesC1]
  n1 <- length(samplesC1)
  E1 <- sum(dist1)/(n1*(n1-1))
  Dis1P2 <- sum(dist1^2)
  E1P2 <- Dis1P2/(n1*(n1-1))
  E1DCross <- (sum(apply(dist1,1,sum)^2)-Dis1P2)/(n1*(n1-1)*(n1-2))
  VarD1 <- E1P2 - E1^2
  CovE1Cross <- E1DCross - E1^2
  VarEta1 <- 4*((n1*(n1-1))/2*VarD1 + n1*(n1-1)*(n1-2)*CovE1Cross)/((n1*(n1-2))^2)
  
  dist2 <- dist[samplesC2,samplesC2]
  n2 <- length(samplesC2)
  E2 <- sum(dist2)/(n2*(n2-1))
  Dis2P2 <- sum(dist2^2)
  E2P2 <- Dis2P2/(n2*(n2-1))
  E2DCross <- (sum(apply(dist2,1,sum)^2)-Dis2P2)/(n2*(n2-1)*(n2-2))
  VarD2 <- E2P2 - E2^2
  CovE2Cross <- E2DCross - E2^2
  VarEta2 <- 4*((n2*(n2-1))/2*VarD2 + n2*(n2-1)*(n2-2)*CovE2Cross)/((n2*(n2-2))^2)
  
  vartotal <- VarEta2 + VarEta1;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  pvalue=2*(1-pnorm(abs(zscore)))
  
  nAll = n1 + n2
  EAll <- sum(dist)/(nAll*(nAll-1))
  DisAllP2 <- sum(dist^2)
  EAllP2 <- DisAllP2/(nAll*(nAll-1))
  EAllDCross <- (sum(apply(dist,1,sum)^2)-DisAllP2)/(nAll*(nAll-1)*(nAll-2))
  CovEAllCross <- EAllDCross - EAll^2
  lambda = n1/(n1+n2)
  
  myvartotal = 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  return(list(E1=E1,E2=E2,
              zscore=zscore,     
              VarEta1=VarEta1,VarEta2=VarEta2, #EVA 1
              sdtotal=sqrt(abs(vartotal)),#EVA 1 sd 
              pvalue=pvalue))
  
}

# function to adjust p values
adjustpval <- function(df) {
  # subset by iteration
  mylist <- split(df, (seq(nrow(df))-1) %/% 50) 
  # apply adjusted pval to each df in list
  for (j in 1:length(mylist)) {
    columnToAdd = p.adjust(mylist[[j]]$pvalue, method='BH')
    mylist[[j]]$pvaluadj = columnToAdd
  }
  return(mylist)
}

# function to generate data with specified zeros and run eva on
beva_path <- function(normal = NULL, nonzero_d, zero_d, nonzero_n, 
                      zero_n, metric) {
  
  exprsdata <- log10(exprsdata)
  if (!is.null(normal)) {
    df <- as.data.frame(exprsdata)
    normal <- df[phenotypes == 0]
    new_normal <- cbind(normal, normal)
    new_phenotypes <- rep(0:1, each = 22)
    dfm <- as.matrix(new_normal)
    rownames(dfm) <- rownames(df)
    newdf <- as.data.frame(dfm)
    
  } else {
    newdf <- as.data.frame(exprsdata)
    new_phenotypes <- phenotypes
  }
  
  d = NULL
  for (i in 1:length(hallmark_pathways)) {
    
    
    pathway = hallmark_pathways[i]
    pathway_name = names(pathway)
    pathway_genes = pathway[[1]]
    df_path <- newdf[(row.names(newdf) %in% pathway_genes), ]
    
    # prune pathway min genes = 5
    numgenes <- as.vector(rownames(df_path))
    if (length(numgenes) <= 5) { 
      next
    }
    
    # add equivalent random zeros to each phenotype by pathway
    disease <- df_path[phenotypes == 1]
    dfdisease <- as.data.frame(lapply(disease, function(cc) cc[sample(c(TRUE, 
                                                                        NA), prob = c(nonzero_d, zero_d), size = length(cc), 
                                                                      replace = TRUE)]))
    dfdisease[is.na(dfdisease)] <- 0
    
    normal <- df_path[phenotypes == 0]
    dfnormal <- as.data.frame(lapply(normal, function(cc) cc[sample(c(TRUE, 
                                                                      NA), prob = c(nonzero_n, zero_n), size = length(cc), 
                                                                    replace = TRUE)]))
    dfnormal[is.na(dfnormal)] <- 0
    path_df <- cbind(dfnormal, dfdisease)
    path_dfm <- as.matrix(path_df)
    rownames(path_dfm) <- rownames(df_path)
    df_path <- as.data.frame(path_dfm)
    
    # distance calculations
    if (metric == "kendall_tau") { 
      distAll <- GSReg.kendall.tau.distance(as.matrix(df_path))
    } else {
      distAll <- distance(t(df_path), method = metric)
    }
    
    colnames(distAll) = toupper(colnames(distAll))
    rownames(distAll) = toupper(rownames(distAll))
    stats <- eva_stats(df = df_path, dist = distAll, pheno = new_phenotypes)
    
    # percent zeros
    logical <- ifelse(newdf == 0, 1,0)
    logical_pathway <- logical[(row.names(logical) %in% pathway_genes),]
    #number of zeros
    zeros <- sum(logical_pathway  == 1) 
    #number of nonzeros
    nonzeros <- sum(logical_pathway  == 0) 
    total <- (zeros + nonzeros)
    percent_zeros <- zeros/total
    percent_nonzeros <- nonzeros/total
    #length of genes in pathway
    genes <- length(pathway_genes)
    
    d = rbind(d, data.frame(pathway_name, stats, percent_zeros, zero_d, zero_n, genes, metric))
  }
  
  return(d)
  
}

################################################# 
### code chunk: run data
################################################# 


# get total number of zeros for each round, report average false positives
# don't need to compare for no signal = all false positives
d = ""
for (i in 1:length(dlist)) {
  distance_metric = dlist[i]
  
  x <- replicate(100, {
    mm <- beva_path(normal = TRUE, 0.6, 0.4, 0.6, 0.4, distance_metric)
    mm$pvaluadj = p.adjust(mm$pvalue, method='BH')
    nrow(mm[mm$pvaluadj < 0.05,])
  })
  
  df <- data.frame(x)
  colnames(df) <- distance_metric
  d <- cbind(d, data.frame(df))
}  

df2 <- data.frame(sapply(d, function(x) as.numeric(as.character(x))))
#write.csv(df2, "EVA_alldist_falsepositives_100_count_40perc.csv")


# data with signal
d = ""
for (i in 34:length(dlist)) {
  distance_metric = dlist[i]
  # true positives - no zeros 
  truth <- beva_path(normal = NULL, 1, 0, 1, 0, distance_metric)
  truth$pvaluadj = p.adjust(truth$pvalue, method='BH')
  sig_truth <- truth[truth$pvaluadj < 0.05,]
  truth_total <- nrow(sig_truth)
  
  x <- replicate(2, {
    
    mm <- beva_path(normal = NULL, 0.8, 0.2, 0.8, 0.2, distance_metric)
    mm$pvaluadj = p.adjust(mm$pvalue, method='BH')
    sig <- mm[mm$pvaluadj < 0.05,]
    # number of true positives
    total <- nrow(sig)
    trues <- nrow(sig[sig$pathway_name %in% sig_truth$pathway_name,])
    falses <- total - trues
    
    rbind(trues, falses)
    })
  
  df <- data.frame(x)
  df <- as.data.frame(t(df))
  df$truth_total <- truth_total
  df$metric <- distance_metric
  d <- rbind(d, data.frame(df))
}

# read in data 

d1 <- read.csv("EVA_alldist_falsepositives_100_count_20perc.csv", header = TRUE, sep = ",")
d1 <- d1[-c(1)]
d2 <- read.csv("EVA_alldist_falsepositives_100_count_40perc.csv", header = TRUE, sep = ",")
d2 <- d2[-c(1)]
d3 <- read.csv("EVA_alldist_falsepositives_100_count_60perc.csv", header = TRUE, sep = ",")
d3 <- d3[-c(1)]
d4 <- read.csv("EVA_alldist_falsepositives_100_count_80perc.csv", header = TRUE, sep = ",")
d4 <- d4[-c(1)]

d1 <- as.data.frame(colSums(d1))
d2 <- as.data.frame(colSums(d2))
d3 <- as.data.frame(colSums(d3))
d4 <- as.data.frame(colSums(d4))

d1$metric <- rownames(d1)
d2$metric <- rownames(d2)
d3$metric <- rownames(d3)
d4$metric <- rownames(d4)

d1$perc <- "20"
d2$perc <- "40"
d3$perc <- "60"
d4$perc <- "80"

common_colnames <- c("total", "metric", "zeros")
colnames(d1) <- common_colnames
colnames(d2) <- common_colnames
colnames(d3) <- common_colnames
colnames(d4) <- common_colnames
dat <- rbind(d1,d2,d3,d4)
dat <- dat[!dat$metric == "d",]

pdf("EVA_alldist_falsepositives_100perm.pdf")
ggplot(dat, aes(x = reorder(metric,-total), y = total, fill = zeros)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  ggtitle("False positives") +
  scale_x_discrete(name="metric") +
  scale_y_discrete(name="Number of false positives") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# read in data - signal

d1 <- read.csv("EVA_alldist_signal_100_count_20perc.csv", header = TRUE, sep = ",")
d1 <- d1[-c(1)]
d2 <- read.csv("EVA_alldist_signal_100_count_40perc.csv", header = TRUE, sep = ",")
d2 <- d2[-c(1)]
d3 <- read.csv("EVA_alldist_signal_100_count_60perc.csv", header = TRUE, sep = ",")
d3 <- d3[-c(1)]
d4 <- read.csv("EVA_alldist_signal_100_count_80perc.csv", header = TRUE, sep = ",")
d4 <- d4[-c(1)]

true1 <- aggregate(trues ~ metric, d1, sum)
true1$zeros <- "0.2"
true2 <- aggregate(trues ~ metric, d2, sum)
true2$zeros <- "0.4"
true3 <- aggregate(trues ~ metric, d3, sum)
true3$zeros <- "0.6"
true4 <- aggregate(trues ~ metric, d4, sum)
true4$zeros <- "0.8"
trues <- rbind(true1,true2,true3,true4)

b <- ggplot(trues, aes(x = reorder(metric,-trues), y = trues, fill = zeros)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  ggtitle("True Positives") +
  scale_x_discrete(name="metric") +
  scale_y_discrete(name="Number of true positives") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

false1 <- aggregate(falses ~ metric, d1, sum)
false1$zeros <- "0.2"
false2 <- aggregate(falses ~ metric, d2, sum)
false2$zeros <- "0.4"
false3 <- aggregate(falses ~ metric, d3, sum)
false3$zeros <- "0.6"
false4 <- aggregate(falses ~ metric, d4, sum)
false4$zeros <- "0.8"
falses <- rbind(false1, false2, false3, false4)

c <- ggplot(falses, aes(x = reorder(metric,-falses), y = falses, fill = zeros)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  ggtitle("False Positives") +
  scale_y_discrete(name="Number of false positives") +
  scale_x_discrete(name="metric") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tot1 <- aggregate(truth_total ~ metric, d1, mean)

a <- ggplot(tot1, aes(x = reorder(metric,-truth_total), y = truth_total, fill=metric)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  ggtitle("Significant pathways w/o zeros") +
  scale_x_discrete(name="metric") +
  scale_y_discrete(name="Number of significant pathways") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=FALSE)

pdf("EVA_alldist_signal_100perm.pdf")
a
b
c
dev.off()
