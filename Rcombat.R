# this is code to replicate the anlayses and figures from my paper
# "Expression variation analysis for tumor heterogeneity in single-cell RNA-sequencing data"
# figure 5A, B, and C
# code developed by Emily Davis

###################################################
### code chunk 1: r init
###################################################
sink(file = "rout.log")
library(ggplot2)
library(mclust)
library(monocle)
library(reshape2)
library(dplyr)
library(GSBenchMark)
library(GSReg)
library(sva)

###################################################
### code chunk 2: load data and pathways
###################################################
# object name cancer
load(file = "magic_cancer_tcga.rda")

# load pathways
load(file = "hallmark_pathways.rda")

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

###################################################
### code chunk : combat and EVA subtypes
###################################################

cdat <- cancer[,pData(cancer)$lymph_node == "NLN"]

# subset for combat
basal <- cdat[,pData(cdat)$subtype == "Basal"]
classical <- cdat[,pData(cdat)$subtype == "Classical"]
atypical <- cdat[,pData(cdat)$subtype == "Atypical"]

basalexprs <- exprs(basal)
classicalexprs <- exprs(classical)
atypicalerexprs <- exprs(atypical)

# combat on basal and classical
pheno = pData(basal)
batch = pheno$patient
modcancer = model.matrix(~as.factor(subtype), data=pheno)
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=basalexprs, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

pheno = pData(classical)
batch = pheno$patient
modcombat = model.matrix(~1, data=pheno)
combat_edata2 = ComBat(dat=classicalexprs, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# rejoin data
combat_dat <- cbind(combat_edata, combat_edata2, atypicalerexprs)

print("combat done")

# EVA
subtype_pathways <- function(pathway, subtype1, subtype2) {
    #E1=0, E2=1
    # subset cancer cells
    pdat <- pData(cdat)
    samples <- pdat[pdat$subtype==subtype1 | pdat$subtype==subtype2,]
    edata <- combat_dat[,rownames(samples)]
    phenoc <- samples
    
    phenoc$compare <- ifelse(phenoc$subtype == subtype1, 0,1)
    
    #eva
    VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=edata, pathways=pathway, phenotypes=as.factor(phenoc$compare))
    
    pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
    kendall <- lapply(VarAnKendallV, function (x) x[1:20])
    dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
    rownames(dd) <- make.unique(names(kendall[[1]]))
    colnames(dd) <- names(kendall)
    dd <- as.data.frame(t(dd))
    dd$pathway_name <- rownames(dd)
    
    dd$metric <- "kendall_tau"
    dd$group1 <- subtype1
    dd$group2 <- subtype2
    return(dd)
}

s1 <- subtype_pathways(pathway = hallmark_pathways, subtype1 = "Atypical", subtype2 = "Basal")
s2 <- subtype_pathways(pathway = nnmf, subtype1 = "Atypical", subtype2 = "Basal")
s3 <- subtype_pathways(pathway = hallmark_pathways, subtype1 = "Classical", subtype2 = "Atypical")
s4 <- subtype_pathways(pathway = nnmf, subtype1 = "Classical", subtype2 = "Atypical")
s5 <- subtype_pathways(pathway = hallmark_pathways, subtype1 = "Basal", subtype2 = "Classical")
s6 <- subtype_pathways(pathway = nnmf, subtype1 = "Basal", subtype2 = "Classical")

all_new <- rbind(s1,s2,s3,s4,s5,s6)
all_new$pvaluadj <- p.adjust(all_new$pvalue,method='BH')
#write.csv(all_new, file = "eva_magic_hnscc_primary_combat_cancer_original_subtypes_results.csv")

common_cols <- c("estat", "pathway", "subtype")
E1 <- all_new[c(1,21,23)]
colnames(E1) <- common_cols
mat <- acast(E1, pathway ~ subtype , value.var='estat')
mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

df = data.frame(subtype = rep(c("Atypical", "Basal", "Classical"), each = 1))

ha1 = HeatmapAnnotation(df = df, col = list(subtype = c('Atypical'="#F8766D",
'Basal'="#7CAE00",
"Classical" = "#00BFC4")),
height = unit(0.5, "cm"))

pdf("tcga_original_subtypes_primary_HNSCC_magic_eva_complexheatmap_scale_new.pdf",width=15,height=10, paper = "USr")
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"), top_annotation = ha1, row_names_side = "left",
row_dend_side = "right", width = unit(60, "mm"), row_names_gp = gpar(fontsize = 8),
column_title = "original TCGA Subtypes")
dev.off()

###################################################
### code chunk : Inter-patient basal
###################################################

basal_pathways <- function(pathway,patient1,patient2) {
    
    pdat <- pData(basal)
    samples <- pdat[pdat$patient==patient1 | pdat$patient==patient2,]
    basalexprs <- exprs(basal)
    edata <- basalexprs[,rownames(samples)]
    phenoc <- samples
    
    phenoc$compare <- ifelse(phenoc$patient == patient1, 0,1)
    
    #eva
    VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=edata, pathways=pathway, phenotypes=as.factor(phenoc$compare))
    
    pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
    kendall <- lapply(VarAnKendallV, function (x) x[1:20])
    dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
    rownames(dd) <- make.unique(names(kendall[[1]]))
    colnames(dd) <- names(kendall)
    dd <- as.data.frame(t(dd))
    dd$pathway_name <- rownames(dd)
    
    dd$metric <- "kendall_tau"
    dd$group1 <- patient1
    dd$group2 <- patient2
    return(dd)
    
    
}

s1 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN18", patient2 = "HN25")
s2 <- basal_pathways(pathway = nnmf, patient1 = "HN18", patient2 = "HN25")
s3 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN16", patient2 = "HN18")
s4 <- basal_pathways(pathway = nnmf, patient1 = "HN16", patient2 = "HN18")
s5 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN18", patient2 = "HN17")
s6 <- basal_pathways(pathway = nnmf, patient1 = "HN18", patient2 = "HN17")
s7 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN18", patient2 = "HN22")
s8 <- basal_pathways(pathway = nnmf, patient1 = "HN18", patient2 = "HN22")
s9 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN18", patient2 = "HN28")
s10 <- basal_pathways(pathway = nnmf, patient1 = "HN18", patient2 = "HN28")
s11 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN18", patient2 = "HN5")
s12 <- basal_pathways(pathway = nnmf, patient1 = "HN18", patient2 = "HN5")
s13 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN25", patient2 = "HN16")
s14 <- basal_pathways(pathway = nnmf, patient1 = "HN25", patient2 = "HN16")
s15 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN25", patient2 = "HN17")
s16 <- basal_pathways(pathway = nnmf, patient1 = "HN25", patient2 = "HN17")
s17 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN25", patient2 = "HN22")
s18 <- basal_pathways(pathway = nnmf, patient1 = "HN25", patient2 = "HN22")
s19 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN25", patient2 = "HN28")
s20 <- basal_pathways(pathway = nnmf, patient1 = "HN25", patient2 = "HN28")
s21 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN25", patient2 = "HN5")
s22 <- basal_pathways(pathway = nnmf, patient1 = "HN25", patient2 = "HN5")
s23 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN16", patient2 = "HN17")
s24 <- basal_pathways(pathway = nnmf, patient1 = "HN16", patient2 = "HN17")
s25 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN16", patient2 = "HN22")
s26 <- basal_pathways(pathway = nnmf, patient1 = "HN16", patient2 = "HN22")
s27 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN16", patient2 = "HN28")
s28 <- basal_pathways(pathway = nnmf, patient1 = "HN16", patient2 = "HN28")
s29 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN16", patient2 = "HN5")
s30 <- basal_pathways(pathway = nnmf, patient1 = "HN16", patient2 = "HN5")
s31 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN17", patient2 = "HN22")
s32 <- basal_pathways(pathway = nnmf, patient1 = "HN17", patient2 = "HN22")
s33 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN17", patient2 = "HN28")
s34 <- basal_pathways(pathway = nnmf, patient1 = "HN17", patient2 = "HN28")
s35 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN17", patient2 = "HN5")
s36 <- basal_pathways(pathway = nnmf, patient1 = "HN17", patient2 = "HN5")
s37 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN22", patient2 = "HN28")
s38 <- basal_pathways(pathway = nnmf, patient1 = "HN22", patient2 = "HN28")
s39 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN22", patient2 = "HN5")
s40 <- basal_pathways(pathway = nnmf, patient1 = "HN22", patient2 = "HN5")
s41 <- basal_pathways(pathway = hallmark_pathways, patient1 = "HN28", patient2 = "HN5")
s42 <- basal_pathways(pathway = nnmf, patient1 = "HN28", patient2 = "HN5")

all_basal <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,
s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42)
all_basal$pvaluadj <- p.adjust(all_basal$pvalue,method='BH')
write.csv(all_basal, file = "eva_magic_hnscc_primary_combat_basal_interpatient_results.csv")

df <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14)
common_cols <- c("estat", "pathway", "patient")
E2 <- df[c(2,21,24)]
colnames(E2) <- common_cols
mat <- acast(E2, pathway ~ patient , value.var='estat')
mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

df = data.frame(subtype = rep(c("Basal"), each = 7))

ha1 = HeatmapAnnotation(df = df, col = list(subtype = c('Basal'="#7CAE00")),
height = unit(0.5, "cm"))

pdf("hallmark_tcga_basal_primary_HNSCC_magic_eva_complexheatmap_scale.pdf",width=15,height=10, paper = "USr")
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"),
top_annotation = ha1,
row_names_side = "left",
row_dend_side = "right", width = unit(60, "mm"), row_names_gp = gpar(fontsize = 8),
column_title = "Basal primaries")
dev.off()
###################################################
### code chunk : Inter-patient fibroblasts from basal patients
###################################################

# EVA on fibroblasts
fibro_pathways <- function(pathway,patient1,patient2) {
    
    fdat <- cds[,pData(cds)$myCell_Type == "NLN-Fibroblast"]
    pdat <- pData(fdat)
    samples <- pdat[pdat$patient==patient1 | pdat$patient==patient2,]
    fexprs <- exprs(fdat)
    edata <- fexprs[,rownames(samples)]
    phenoc <- samples
    phenoc$compare <- ifelse(phenoc$patient == patient1, 0,1)
    
    VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=edata, pathways=pathway, phenotypes=as.factor(phenoc$compare))
    
    pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
    kendall <- lapply(VarAnKendallV, function (x) x[1:20])
    dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
    rownames(dd) <- make.unique(names(kendall[[1]]))
    colnames(dd) <- names(kendall)
    dd <- as.data.frame(t(dd))
    dd$pathway_name <- rownames(dd)
    
    dd$metric <- "kendall_tau"
    dd$group1 <- patient1
    dd$group2 <- patient2
    return(dd)
}

s1 <- fibro_pathways(pathway = myeloid, patient1 = "HN18", patient2 = "HN25")
s2 <- fibro_pathways(pathway = myeloid, patient1 = "HN18", patient2 = "HN16")
s3 <- fibro_pathways(pathway = myeloid, patient1 = "HN18", patient2 = "HN17")
s4 <- fibro_pathways(pathway = myeloid, patient1 = "HN18", patient2 = "HN22")
s5 <- fibro_pathways(pathway = myeloid, patient1 = "HN18", patient2 = "HN28")
s6 <- fibro_pathways(pathway = myeloid, patient1 = "HN5", patient2 = "HN18")
s7 <- fibro_pathways(pathway = myeloid, patient1 = "HN25", patient2 = "HN16")
s8 <- fibro_pathways(pathway = myeloid, patient1 = "HN25", patient2 = "HN17")
s9 <- fibro_pathways(pathway = myeloid, patient1 = "HN25", patient2 = "HN22")
s10 <- fibro_pathways(pathway = myeloid, patient1 = "HN25", patient2 = "HN28")
s11 <- fibro_pathways(pathway = myeloid, patient1 = "HN25", patient2 = "HN5")
s12 <- fibro_pathways(pathway = myeloid, patient1 = "HN16", patient2 = "HN17")
s13 <- fibro_pathways(pathway = myeloid, patient1 = "HN16", patient2 = "HN22")
s14 <- fibro_pathways(pathway = myeloid, patient1 = "HN16", patient2 = "HN28")
s15 <- fibro_pathways(pathway = myeloid, patient1 = "HN16", patient2 = "HN5")
s16 <- fibro_pathways(pathway = myeloid, patient1 = "HN17", patient2 = "HN22")
s17 <- fibro_pathways(pathway = myeloid, patient1 = "HN17", patient2 = "HN28")
s18 <- fibro_pathways(pathway = myeloid, patient1 = "HN17", patient2 = "HN5")
s19 <- fibro_pathways(pathway = myeloid, patient1 = "HN22", patient2 = "HN28")
s20 <- fibro_pathways(pathway = myeloid, patient1 = "HN22", patient2 = "HN5")
s21 <- fibro_pathways(pathway = myeloid, patient1 = "HN28", patient2 = "HN5")

all_fibro <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21)
all_fibro$pvaluadj <- p.adjust(all_fibro$pvalue,method='BH')
write.csv(all_fibro, file = "eva_magic_hnscc_basal_fibroblast_interpatient_results.csv")

df <- rbind(s1,s2,s3,s4,s5,s6,s11)
common_cols <- c("estat", "pathway", "patient")
E2 <- df[c(2,21,24)]
colnames(E2) <- common_cols

mat <- acast(E2, pathway ~ patient , value.var='estat')

mats <- t(apply(mat,1,scale))
colnames(mats) <- colnames(mat)

df = data.frame(fibroblast = c("low", "low", "low", "low" ,"high" ,"high", "low"))

ha1 = HeatmapAnnotation(df = df, col = list(fibroblast = c('high'="red",
'low'="blue")),
height = unit(0.5, "cm"))

pdf("tcga_basal_fibroblasts_HNSCC_magic_eva_complexheatmap_scale.pdf",width=15,height=10, paper = "USr")
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"),
top_annotation = ha1,
row_names_side = "left",
row_dend_side = "right", width = unit(60, "mm"), row_names_gp = gpar(fontsize = 8),
column_title = "Basal patient fibroblast - Myeloid")
dev.off()
