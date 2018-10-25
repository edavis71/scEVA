# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# convert nanostring file to list of pathways
# code developed by Emily Davis

###################################################
### code chunk: pathways
###################################################

# read pathway info from nanostring
# be sure to remove header images prior to loading
pan <- read.csv('mmc4.csv',header=T, sep=",")

# list of annotation categories
cat <- as.vector(colnames(pan))
rownames(pan) <- pan$Gene

# tcell
mylist <- NULL
for (i in 1:length(cat)) {
  trues <- pan[cat[i]]
  trues <- trues[ rowSums(trues != "") != 0,]
  trues <- as.character(trues)
  name <- cat[i]
  print(name)
  mylist[[name]] <- trues
}

# nanostring
mylist <- NULL
for (i in 1:length(cat)) {
  trues <- pan[which(pan[cat[i]] == "+"),]
  genes <- toupper(rownames(trues))
  name <- cat[i]
  print(name)
  mylist[[name]] <- genes
}

tcell_activation <- mylist
save(tcell_activation, file = "tcell_activation_pathways.rda")

