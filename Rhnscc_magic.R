# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 5C and D
# code developed by Emily Davis

###################################################
### code chunk 1: r init
###################################################
library(ggplot2)
library(mclust)
library(monocle)
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(GSBenchMark)
library(GSReg)

###################################################
### code chunk 2: load data
###################################################

load("magic_hnscc_cds.rda")

# load pathways
load(file = "hallmark_pathways.rda")

emt = list(EEMT = c("CLDN4", "DYNC1LI2", "IRF6", "CTNND1", "HOOK1", "GALNT3", "MYO5B", "PRSS8", "GPR56", "CNOT1", "GRHL2", "ERBB3", "ATP8B1", "F11R", "OCLN", "CDS1", "MAP7", "CGN", "SPINT1", "ESRP1", "AP1G1", "MARVELD3", "MARVELD2", "ESRP2", "CDH1"), MEMT = c("SPOCK1", "SYT11", "NAP1L3", "CALD1", "BNC2", "ANTXR1", "DACT1", "GPC6", "CMTM3", "FSTL1", "MSRB3", "GYPC", "HTRA1", "AXL", "ANGPTL2", "ZEB2", "PMP22", "COL6A1", "LOXL2", "PCOLCE", "PDGFRB", "EMP3", "COL6A2", "MMP2", "CNRIP1", "POSTN", "ADAM12", "COL3A1", "COL10A1", "SULF1", "NID2", "LRRC15", "COL8A1", "FAP", "VCAN", "ITGA11", "INHBA", "SPARC", "OLFML2B", "AEBP1", "ADAMTS12", "ADAMTS2", "COL1A1", "COL1A2", "COL5A1", "THBS2", "COL5A2", "COL6A3", "FBN1", "CDH2", "FN1", "VIM"))

nnmf = list(
CELL_CYCLE = c("TK1","HMGB2","ZWINT","MAD2L1","TUBA1B","STMN1","KIF22","CKS1B","H2AFZ","CENPW","CDC20","DTYMK","UBE2C","UBE2T","NUSAP1","RRM2","BIRC5","RNASEH2A","PCNA","TUBB","KPNA2","ASF1B","TRIP13","CCNB1","TPX2","CCNB2","TYMS","PTTG1","KIAA0101","GMNN","DNAJC9","CCNA2","CKS2","MLF1IP","VRK1","CENPM","PRC1","SPAG5","TOP2A","AURKB","FEN1","TMEM106C","RRM1","RFC4","MCM7","CDKN3","NUDT1","PBK","MELK","ANLN","MCM5","PLK1","GGH","MCM4","CENPN","TMPO","CDCA3","DEK","RPA2","KIF2C","CDK1","CDCA5","LSM4","KNSTRN","TUBG1","SMC4","CSE1L","UHRF1","RANBP1","CDCA8","MCM2","RFC2","HMGN2","ATAD2","HAT1","PKMYT1","SIVA1","FANCI","ECT2","POLE3","WDR34","MCM3","NCAPG2","TUBB6","NCAPD2","GINS2","TIMELESS","RAD51","CMC2","OIP5","TUBB4B","APOBEC3B","ORC6","C19orf48","SNRNP25","RFC3","TROAP","EBP","DKC1","H2AFV"),
PEMT = c("SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC","MMP2","EMP3","INHBA","LAMB3","VIM","SEMA3C","PRKCDBP","ANXA5","DHRS7","ITGB1","ACTN1","CXCR7","ITGB6","IGFBP7","THBS1","PTHLH","TNFRSF6B","PDLIM7","CAV1","DKK3","COL17A1","LTBP1","COL5A2","COL1A1","FHL2","TIMP3","PLAU","LGALS1","PSMD2","CD63","HERPUD1","TPM1","SLC39A14","C1S","MMP1","EXT2","COL4A2","PRSS23","SLC7A8","SLC31A2","ARPC1B","APP","MFAP2","MPZL1","DFNA5","MT2A","MAGED2","ITGA6","FSTL1","TNFRSF12A","IL32","COPB2","PTK7","OCIAD2","TAX1BP3","SEC13","SERPINH1","TPM4","MYH9","ANXA8L1","PLOD2","GALNT2","LEPREL1","MAGED1","SLC38A5","FSTL3","CD99","F3","PSAP","NMRK1","FKBP9","DSG2","ECM1","HTRA1","SERINC1","CALU","TPST1","PLOD3","IGFBP3","FRMD6","CXCL14","SERPINE2","RABAC1","TMED9","NAGK","BMP1","ESYT1","STON2","TAGLN","GJA1"),
EPI_DIFF1 = c("IL1RN","SLPI","CLDN4","S100A9","SPRR1B","PVRL4","RHCG","SDCBP2","S100A8","APOBEC3A","GRHL1","SULT2B1","ELF3","KRT16","PRSS8","MXD1","S100A7","KRT6B","LYPD3","TACSTD2","CDKN1A","KLK11","GPRC5A","KLK10","TMBIM1","PLAUR","CLDN7","DUOXA1","PDZK1IP1","NCCRP1","IDS","PPL","ZNF750","EMP1","CLDN1","CRB3","CYB5R1","DSC2","S100P","GRHL3","SPINT1","SDR16C5","SPRR1A","WBP2","GRB7","KLK7","TMEM79","SBSN","PIM1","CLIC3","MALAT1","TRIP10","CAST","TMPRSS4","TOM1","A2ML1","MBOAT2","LGALS3","ERO1L","EHF","LCN2","YPEL5","ALDH3B2","DMKN","PIK3IP1","CEACAM6","OVOL1","TMPRSS11E","CD55","KLK6","SPRR2D","NDRG2","CD24","HIST1H1C","LY6D","CLIP1","HIST1H2AC","BNIPL","QSOX1","ECM1","DHRS3","PPP1R15A","TRIM16","AQP3","IRF6","CSTA","RAB25","HOPX","GIPC1","RAB11FIP1","CSTB","KRT6C","PKP1","JUP","MAFF","DSG3","AKTIP","KLF3","HSPB8","H1F0"),
EPI_DIFF2 =
c("LY6D","KRT16","KRT6B","LYPD3","KRT6C","TYMP","FABP5","SCO2","FGFBP1","JUP","IMP4","DSC2","TMBIM1","KRT14","C1QBP","SFN","S100A14","RAB38","GJB5","MRPL14","TRIM29","ANXA8L2","KRT6A","PDHB","AKR1B10","LAD1","DSG3","MRPL21","NDUFS7","PSMD6","AHCY","GBP2","TXN2","PSMD13","NOP16","EIF4EBP1","MRPL12","HSD17B10","LGALS7B","THBD","EXOSC4","APRT","ANXA8L1","ATP5G1","S100A2","TBRG4","MAL2","NHP2L1","DDX39A","ZNF750","UBE2L6","WDR74","PPIF","PRMT5","VSNL1","VPS25","SNRNP40","ADRM1","NDUFS8","TUBA1C","TMEM79","UQCRFS1","EIF3K","NME2","PKP3","SERPINB1","RPL26L1","EIF6","DSP","PHLDA2","S100A16","LGALS7","MT1X","UQCRC2","EIF3I","MRPL24","CCT7","RHOV","ECE2","SSBP1","POLDIP2","FIS1","CKMT1A","GJB3","NME1","MRPS12","GPS1","ALG3","MRPL20","EMC6","SRD5A1","PA2G4","ECSIT","MRPL23","NAA20","HMOX2","COA4","DCXR","PSMD8","WBSCR22"),
STRESS =
c("FOS","ATF3","NR4A1","DUSP1","ZFP36","PPP1R15A","SGK1","EGR1","ZC3H12A","JUNB","FOSB","IER2","NFKBIA","NFKBIZ","HBEGF","BTG2","SOD2","CDKN1A","NCOA7","JUN","MYC","SERTAD1","CCNL1","RND3","PLK2","SOCS3","DNAJB1","DUSP2","TSC22D1","KLF10","GADD45B","PMAIP1","MAFF","ERRFI1","SLC38A2","IRF1","TOB1","ID2","KLF6","DNAJA1","TNFAIP3","BHLHE40","NXF1","FOSL1","IER3","DUSP6","HCAR2","IL8","CYR61","EFNA1","C1R","PHLDA2","DNAJB14","MCL1","HERPUD1","ADRB2","EIF4A3","TACSTD2","ID1","ETS2","CD74","TRIB1","SLC20A1","LOC284454","EIF1","CXCL2","BRD2","RASD1","LDLR","EGR2","TFRC","ADM","TGIF1","HLA-DRB1","OSR2","SAA1","ELF3","CLK1","PER2","KLF4","GPNMB","MXD1","UBC","HLA-DRA","SLC3A2","OVOL1","HIST1H2BK","DDX3X","LAMB3","ZNF622","TUBB2A","ZFAND5","IRF6","TNF","BTG1","LMNA","MAP1LC3B","TSC22D3","PLK3","KLHL21"),
HYPOXIA = c("NDRG1","IGFBP3","PTHLH","EGLN3","BNIP3","NDUFA4L2","ERO1L","P4HA1","SLC2A1","ENO2","HK2","PGF","LDHA","PGK1","PDK1","DHRS3","DDIT4","PVRL4","GPNMB","BIK","GJB6","C4orf3","IGFBP2","FAM162A","GPI","LPIN3","PLAU","ADM","ANGPTL4","DARS","NUPR1","SERPINE1","PGAM1","ALDOA","DAAM1","CXADR","SEMA4B","CA9","CIB1","SPRR1B","PLIN2","WSB1","HILPDA","NOL3","PFKFB3","IFNGR1","H1F0","KDM3A","BCL6","BNIP3L","ZFP36L1","HLA-E","PIK3IP1","CLK3","POLR1D","BTG1","NPC2","LAMP2","DSG2","SAT1","AK4","SMS","FRMD6","CLDND1","ACP5","AP1G2","TPI1","PLAUR","BCL10","TMEM59","HAS3","SERINC1","C1orf43","ENO1","CSDA","PFKP","KLHL24","HIST1H1C","RBPJ","BHLHE40","GAPDH","UPK3BL","LTBP1","P4HA2","HBP1","GRHL1","DDIT3","ANXA1","ITGA5","LOC100862671","PLS3","TSC22D2","GLTP","PLOD2","PERP","MALL","CTNND1","KDM5B","AHNAK2","PNRC1"))

print("pathways loaded")

table(pData(cds)$myCell_Type, pData(cds)$patient)

###################################################
### code chunk 3: Primary NLN vs Met LN all patients EVA
###################################################

do_pathways <- function(pathway, sample) {
  # subset cds by patient ID
  test<- cds[,(pData(cds)$patient==sample)]
  # subset cancer cells
  testc <- test[,(pData(test)$cell_type=="Cancer")]
  #ln =1; nln =0 in pData(testc)[,3]
  # get expression matrix for pt cancer cells
  exprsdata <- as.matrix(exprs(testc))

  #classic eva
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway, phenotypes=as.factor(pData(testc)[,3])) 

  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)

  dd$metric <- "kendall_tau"
  dd$patient <- sample
  dd$group1 <- "NLN"
  dd$group2 <- "LN"
  return(dd)
}
         
a1 <- do_pathways(pathway = hallmark_pathways, sample = "HN5")
a2 <- do_pathways(pathway = hallmark_pathways, sample = "HN20")
a3 <- do_pathways(pathway = hallmark_pathways, sample = "HN25")
a4 <- do_pathways(pathway = hallmark_pathways, sample = "HN26")
a5 <- do_pathways(pathway = hallmark_pathways, sample = "HN28")

b1 <- do_pathways(pathway = emt, sample = "HN5")
b2 <- do_pathways(pathway = emt, sample = "HN20")
b3 <- do_pathways(pathway = emt, sample = "HN25")
b4 <- do_pathways(pathway = emt, sample = "HN26")
b5 <- do_pathways(pathway = emt, sample = "HN28")

c1 <- do_pathways(pathway = nnmf, sample = "HN5")
c2 <- do_pathways(pathway = nnmf, sample = "HN20")
c3 <- do_pathways(pathway = nnmf, sample = "HN25")
c4 <- do_pathways(pathway = nnmf, sample = "HN26")
c5 <- do_pathways(pathway = nnmf, sample = "HN28")

all <- rbind(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,c1,c2,c3,c4,c5)
all$pvaluadj <- p.adjust(all$pvalue,method='BH')
#write.csv(all, file = "eva_magic_hnscc_lnvsnln_all_results.csv")

all$logpvaluadj <- log10(all$pvaluadj)

sig <- all[which(all$pvaluadj < 0.05),]

df <- read.csv('eva_magic_hnscc_lnvsnln_all_results.csv',header=T, sep=",")
rownames(df) <- df$X
df <- df[-c(1)]

common_cols <- c("estat", "pathway", "patient", "celltype")
E1 <- df[c(1,21,23,24)]
colnames(E1) <- common_cols
E2 <- df[c(2,21,23,25)]
colnames(E2) <- common_cols
dat <- rbind(E1,E2)


pdf("intrapatient_NLNvsLN_HNSCC_magic_eva_sig.pdf",width=15,height=10, paper = "USr")
ggplot(dat, aes(celltype, pathway)) +
  geom_tile(aes(fill = estat), color = "white") +
  scale_fill_viridis("Estat", direction = 1) + 
  ggtitle("Significant estats LN vs NLN") + 
  theme(axis.text.x=element_text(angle=90, hjust =1), plot.title = element_text(hjust = 0.5)) +facet_wrap(~patient)
dev.off()

# ggplot heatmap
library(ggdendro)
library(gridExtra)
library(reshape2)

sig <- all[which(all$pvaluadj < 0.05),]

sigpathways <- unique(all$pathway_name)

# ggplot clustered heatmap
heatdat <- acast(all, patient~pathway_name, value.var="logpvaluadj")
data <- apply(heatdat,2,scale)
rownames(data) <- rownames(heatdat)

datat <- t(data)
ord <- as.dendrogram(hclust( dist(datat, method = "euclidean"), method = "ward.D" ))

# Create dendro
dendro.plot <- ggdendrogram(data = ord, rotate = TRUE)

# Preview the plot
print(dendro.plot)

pd <- as.data.frame( data )
pd$patient <- rownames(data)
pd.m <- melt( pd, id.vars = "patient", variable.name = "pathway" )

otter.order <- order.dendrogram(ord)
pd.m$pathway <- factor(x = pd.m$pathway,
                      levels = colnames(data)[otter.order], 
                      ordered = TRUE)


heatmap.plot <- ggplot( pd.m, aes(patient, pathway) ) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis("log10(pvalue)", direction = 1) +
  labs(title = "Intrapatient NLN vs LN HNSCC") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 6), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.position = "top")

grid.newpage()
#pdf("intrapatient_NLNvsLN_HNSCC_magic_eva_heatmap.pdf",width=15,height=10, paper = "USr")
print(heatmap.plot, vp = viewport(x = 0.3, y = 0.5, width = 0.5, height = 1))
print(dendro.plot, vp = viewport(x = 0.8, y = 0.45, width = 0.5, height = 0.93))
#dev.off()


#all samples 
dat$type <- NA
dat$type[(dat$celltype == "NLN")] <- "primary"
dat$type[(dat$celltype == "LN")] <- "met" 
dat$col <- paste(dat$patient, dat$type, sep = "_")

mat <- acast(dat, pathway ~ col , value.var='estat')

mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

df = data.frame(patient = rep(c("HN5", "HN20", "HN25", "HN26", "HN28"), times = 2),
                celltype = rep(c("met", "primary"), each = 5))

rownames(df) <- paste(df$patient, df$celltype, sep = "_")
#set order by colnames in mat
df <- df[sort(colnames(mat)),]

ha1 = HeatmapAnnotation(df = df, col = list(patient = c('HN5'='red',
                                                        'HN20'='orange',
                                                        'HN25'='yellow',
                                                        'HN26'='green',
                                                        'HN28'='blue'),
                                            celltype = c("met" = "dodgerblue4", "primary" = "red")))
                                            
      
pdf("all_intrapatient_NLNvsLN_HNSCC_magic_eva_complexheatmap_scale.pdf",width=15,height=10, paper = "USr")
Heatmap(mats, c("steelblue3", "khaki1", "red1"), top_annotation = ha1, row_names_side = "left", 
        row_dend_side = "right", width = unit(60, "mm"), row_names_gp = gpar(fontsize = 8))
dev.off()

# samples 25 and 26

sub = dat[ (dat$patient == "HN25" | dat$patient == "HN26") ,]

sub <- sub[!sub$pathway %in% c("EEMT", "MEMT"),]
mat <- acast(sub, pathway ~ col , value.var='estat')

mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

df = data.frame(patient = rep(c("HN25", "HN26"), times = 2),
                celltype = rep(c("metastasis", "primary"), each = 2))

rownames(df) <- paste(df$patient, df$celltype, sep = "_")
#set order by colnames in mat
df <- df[sort(colnames(mat)),]

ha1 = HeatmapAnnotation(df = df, col = list(patient = c('HN25'="#7CAE00",
                                                        'HN26'="#C77CFF"),
                                            celltype = c("metastasis" = "#F8766D", "primary" = "#00BFC4")), 
                        height = unit(0.5, "cm"))


pdf("sig_intrapatient_NLNvsLN_HNSCC_magic_eva_complexheatmap_scale_ggcolors.pdf",width=12,height=20, paper = "USr")
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"), top_annotation = ha1, row_names_side = "left", 
        row_dend_side = "right", width = unit(30, "mm"), row_names_gp = gpar(fontsize = 8))
dev.off()


###################################################
### code chunk : tsne HN25 and HN26 cancer cells
###################################################

load("hnscc_cds_nozeros.rda")


dat <- cds[,pData(cds)$cell_type == "Cancer"]
dat <- dat[,pData(dat)$patient == "HN25" | pData(dat)$patient == "HN26" ]

################################
# model prep
dat <- estimateSizeFactors(dat)
dat <- estimateDispersions(dat)

##############################
# cluster cells
disp_table <- dispersionTable(dat)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.001 & dispersion_empirical >= 2 * dispersion_fit)

dat <- setOrderingFilter(dat, unsup_clustering_genes$gene_id)
plot_ordering_genes(dat)

############################
### pca 
plot_pc_variance_explained(dat, return_all = F) # norm_method = 'log',

#Plot PCA with fewer total components
plot_pc_variance_explained(dat, return_all = F, max_components = 25) # norm_method = 'log',
num_dims<-10

###########################
### reduce dimensions
dat <- reduceDimension(dat, max_components = 2, num_dim = num_dims,
                         reduction_method = 'tSNE', verbose = T,reducedModelFormulaStr="~num_genes_expressed")

dat <- clusterCells(dat)
plot_cell_clusters(dat, 1, 2, color="Cluster",cell_size=.1)

plot_cell_clusters(dat, color = "patient", cell_size=1.2)
plot_cell_clusters(dat, color = "myCell_Type", cell_size=1.2)

