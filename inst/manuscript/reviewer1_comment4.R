library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(enrichR)

enrichmentComparison <- function(X,Y){
  eX <- do.call(rbind.data.frame, enrichr(X, databases = c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")))
  eX <- eX$Term[eX$Adjusted.P.value < 0.05]
  eY <- do.call(rbind.data.frame, enrichr(Y, databases = c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")))
  eY <- eY$Term[eY$Adjusted.P.value < 0.05]
  TP = sum(eX %in% eY)
  FP = sum(!eY %in% eX)
  FN = sum(!eX %in% eY)
  c(Shared=TP, DROnly = FP, DEONly=FN)
}

# TREM2
CM <- read.csv('TREM2/Data/GSE130626_umi_counts.csv.gz')
MD <- read.csv('TREM2/Data/GSE130626_cell_info.csv.gz')
GD <- read.csv('TREM2/Data/GSE130626_gene_info.csv.gz')

CM <- CM[!is.na(GD$symbol),]
GD <- GD[!is.na(GD$symbol),]

rownames(CM) <- make.unique(GD$symbol)

MD <- MD[((MD$treatment %in% 'cuprizone') & (MD$trem2_genotype %in% c('WT', 'KO'))),]
CM <- CM[,MD$cell_id]

D <- CreateSeuratObject(CM)
D <- NormalizeData(D)
D <- ScaleData(D)

DE <- FindMarkers(D, ident.1 = 'WT', ident.2 = 'KO', test.use = 'MAST')
DR <- read.csv('TREM2/Results/drTREM2.csv', row.names = 1)

genesDE <- rownames(DE)[((DE$p_val_adj < 0.05) & (abs(DE$avg_log2FC) > 1))]
# Cd52, Lpl, Spp1, Trem2, Ttr
genesDR <- DR$gene[DR$p.adj < 0.05]
# Adssl1, Akr1a1, Aldh2, Aldoa, Anp32b, Anxa2, Anxa5, Aplp2, Apoe, Arpc1b, Ass1, Atox1, Atp5e, Atp5g1, Atp5h, Atp5j, Atp5j2, Atp5k, Atp5l, Atp5mpl, Atp6v1f, Atpif1, Axl, Bola2, Calm1, Capg, Ccl5, Ccl6, Cd52, Cd63, Cd68, Cd72, Cd9, Ch25h, Chchd10, Cops9, Cotl1, Cox5b, Cox6a1, Cox6a2, Cox6b1, Cox6c, Cox7a2, Cox7c, Cox8a, Crip1, Csf2ra, Cst7, Cstb, Ctsa, Ctsb, Ctsd, Ctsl, Ctsz, Cxcl14, Cxcl16, Cybb, Cycs, Dbi, Eef1g, Eif5a, Elob, Fabp3, Fabp5, Fth1, Fxyd5, Gabarap, Gapdh, Gla, Glipr1, Gm10053, Gm10076, Gm2000, Gm49339, Gm8730, Gnas, Gng5, Gpx4, H2-D1, Hint1, Ifi27l2a, Ifitm3, Igf1, Itgax, Lgals3, Lpl, Lrpap1, Lyz2, Mif, Mmp12, Ms4a6c, Myl6, Ndufa1, Ndufa13, Ndufa2, Ndufb2, Npc2, Ost4, Pgam1, Pgk1, Pkm, Pld3, Plin2, Prdx1, Prdx5, Rab7b, Rap2b, Rbx1, Rhoc, Romo1, S100a1, Selenow, Sem1, Spp1, Syngr1, Tceal9, Tmem256, Tmsb10, Tomm7, Tpi1, Tpm4, Trem2, Tyrobp, Uqcr11, Uqcrb, Uqcrh, Vat1, Ybx1

length(intersect(genesDE, genesDR))/length(genesDE)
# % Overlap between DR and DE
# 0.8

enrichmentComparison(genesDE, genesDR)
# Shared DROnly DEONly
#   23    175    122

FC <- DE$avg_log2FC
names(FC) <- rownames(DE)

t.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR], alternative = 'greater')
# Welch Two Sample t-test
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# t = 3.3214, df = 39.453, p-value = 0.0009688
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.2143344       Inf
# sample estimates:
#   mean of x mean of y
# 0.5386503 0.1037793

wilcox.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR])
# Wilcoxon rank sum test with continuity correction
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# W = 1687, p-value = 0.0008137
# alternative hypothesis: true location shift is not equal to 0

plotData <- data.frame(G =names(FC), FC = FC, DR = factor(ifelse(names(FC) %in% genesDR, 'Yes', 'No'), levels = c('Yes', 'No')))
png('TREM2/Results/drde_Trem2.png', width = 800, height = 1000, res = 300)
ggplot(plotData, aes(DR, FC)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  geom_jitter(alpha = .5, pch = 16) +
  xlab('scTenifoldKnk\nDifferentially Regulated') +
  ylab(parse(text = 'log[2]~(Fold-Change)~by~MAST')) +
  labs(title = 'Trem2') +
  theme(plot.title = element_text(face = 2)) +
  geom_abline(slope = 0, intercept = c(-1,1), lty = 2, col = 'red')
dev.off()

# NKX2-1
WT <- readMM('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_matrix.mtx.gz')
KO <- readMM('NKX2-1/Data/GSM3716704_Nkx2-1_mutant_scRNAseq_matrix.mtx.gz')
rownames(WT) <- read.csv('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE)[,2]
rownames(KO) <- read.csv('NKX2-1/Data/GSM3716704_Nkx2-1_mutant_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE)[,2]
colnames(WT) <- readLines('NKX2-1/Data/GSM3716703_Nkx2-1_control_scRNAseq_barcodes.tsv.gz')
colnames(KO) <- readLines('NKX2-1/Data/GSM3716704_Nkx2-1_mutant_scRNAseq_barcodes.tsv.gz')

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')
NKX21 <- merge(WT, KO)

NKX21 <- NormalizeData(NKX21)
NKX21 <- ScaleData(NKX21)

DE <- FindMarkers(NKX21, ident.1 = 'WT', ident.2 = 'KO', test.use = 'MAST')
DR <- read.csv('NKX2-1/Results/dr_GSM3716703.csv', row.names = 1)

genesDE <- rownames(DE)[(abs(DE$avg_log2FC) > 1 & (DE$p_val_adj < 0.05))]
# 2200002D01Rik, Ager, Apoc1, Apoe, Aqp5, Arl4c, Bex2, Cbr2, Ccnd2, Cd24a, Cd2ap, Cd74, Chil1, Clu, Cxcl15, Cym, Cystm1, Foxq1, Fxyd3, Gm42418, Gsta4, Gsto1, H2-Aa, H2-Ab1, Hist1h1e, Igfbp2, Jun, Krt18, Krt19, Krt8, Lamp3, Lgals3, Lgi3, Lpcat1, Lurap1l, Lyz1, Lyz2, Mal, Msln, Napsa, Nkx2-1, Pdzk1ip1, Ppp1r14c, Retnla, Rnase4, S100a11, S100a6, S100g, Sfta2, Sftpa1, Sftpb, Sftpc, Slc34a2, Spp1, Sprr1a, Sprr2a3, Tff1, Tff2, Tgfbi, Thbs1, Tmsb4x, Tnfrsf12a, Trf
genesDR <- DR$gene[DR$p.adj < 0.05]
# 1810011O10Rik, Abca3, Ace, Acot1, Acvrl1, Ager, Alcam, Ank3, Aqp5, Arhgef15, Atp11a, Atp1b1, Atp6v1c2, Atp8a1, AU021092, Avpi1, BC028528, Bcam, Bex2, Bex4, Calcrl, Car8, Cat, Cbr2, Cd74, Cd93, Cdh1, Cdh5, Cebpa, Ces1d, Chchd10, Chil1, Cldn18, Cldn3, Cldn5, Cldn7, Clec14a, Clec1a, Clec2d, Clic3, Crb3, Crlf1, Ctla2a, Ctsc, Ctsh, Cxcl15, Cxx1a, Cystm1, Cyyr1, Dcxr, Dram1, Ecscr, Egfl6, Egfl7, Emb, Emcn, Epas1, Epcam, Etv5, Exosc7, Ezr, Fabp5, Fendrr, Fgfr2, Flt1, Fmo1, Foxf1, Gata2, Gde1, Gimap6, Gng11, Gpihbp1, Gprc5a, H2afj, Hc, Hilpda, Hlx, Hmcn1, Hoxa5, Hpgd, Icam1, Icam2, Impdh1, Irx1, Irx2, Irx3, Klf2, Krt18, Krt8, Lamp3, Ldb2, Lgi3, Lmo2, Lmo7, Lpcat1, Marcks, Matn4, Mbip, Mcam, Mest, Mettl7a1, Mfap2, Mgst1, Muc1, Myct1, Myh14, Myzap, Napsa, Neat1, Nid1, Nkx2-1, Nostrin, Npc2, Npw, Nrarp, Ociad2, Pcdh1, Pcdh17, Pecam1, Pi4k2b, Pla2g1b, Plpp3, Plvap, Ppp1r14a, Ppp1r14c, Prex2, Prss8, Ptprf, Rab25, Ramp2, Rasip1, Rbpjl, Rhoj, Rnase4, S100g, S1pr1, Scn3b, Scn7a, Sdc1, Sdc4, Sema3c, Sept4, Serpinb6b, Sfta2, Sftpa1, Sftpb, Sftpc, Sftpd, Slc1a5, Slc34a2, Slc39a8, Soat1, Sox17, Sparcl1, Spint2, Srgn, Stmn2, Tbx2, Tcf4, Thbd, Tie1, Tmem100, Tmem204, Tmem238, Tmem243, Trf, Tspan15, Tspan7, Vamp8, Wfdc2, Zeb2

length(intersect(genesDE, genesDR))/length(genesDE)
# % Overlap between DR and DE
# 0.3809524

enrichmentComparison(genesDE, genesDR)
# Shared DROnly DEONly
#   17     71     66

FC <- DE$avg_log2FC
names(FC) <- rownames(DE)

t.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR], alternative = 'greater')
# Welch Two Sample t-test
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# t = 7.6433, df = 106.67, p-value = 4.875e-12
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.5274326       Inf
# sample estimates:
#   mean of x  mean of y
# 0.4680584 -0.2056227

wilcox.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR])
# Wilcoxon rank sum test with continuity correction
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# W = 81026, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

plotData <- data.frame(G =names(FC), FC = FC, DR = factor(ifelse(names(FC) %in% genesDR, 'Yes', 'No'), levels = c('Yes', 'No')))
png('NKX2-1/Results/drde_Nkx2-1.png', width = 800, height = 1000, res = 300)
ggplot(plotData, aes(DR, FC)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  geom_jitter(alpha = .5, pch = 16) +
  xlab('scTenifoldKnk\nDifferentially Regulated') +
  ylab(parse(text = 'log[2]~(Fold-Change)~by~MAST')) +
  labs(title = 'Nkx2-1') +
  theme(plot.title = element_text(face = 2)) +
  geom_abline(slope = 0, intercept = c(-1,1), lty = 2, col = 'red')
dev.off()

# HNF4AG
WT <- Read10X_h5('HNF4A-HNF4G/Data/GSM3477499_WT_ScRNAseq_filtered_gene_bc_matrices.h5')
KO <- Read10X_h5('HNF4A-HNF4G/Data/GSM3477500_Hnf4alphagammaDKO_ScRNAseq_filtered_gene_bc_matrices.h5')

WT <- CreateSeuratObject(WT, project = 'WT')
KO <- CreateSeuratObject(KO, project = 'KO')

HNF4AG <- merge(WT, KO)
HNF4AG <- NormalizeData(HNF4AG)
HNF4AG <- ScaleData(HNF4AG)

DE <- FindMarkers(HNF4AG, ident.1 = 'WT', ident.2 = 'KO', test.use = 'MAST')
DR <- read.csv('HNF4A-HNF4G/Results/drHNF4AG.csv', row.names = 1)

genesDR <- DR$gene[DR$p.adj < 0.05]
# Ace2, Anpep, Apoa1, Apoa4, Apob, Apoc3, Atf3, Atp5e, Atp5j2, Atp5k, Atp5l, Atpif1, Cd36, Cd74, Cox17, Cox5b, Cyp4v3, Dgat1, Dnase1, Eef1a1, Eef2, Fau, Ggt1, Gng5, Gpx1, Guca2b, H2-Aa, H2-Ab1, H2-Q1, H2-Q2, Hnf4a, Hnf4g, Hsp90ab1, Hspa8, Itm2b, Krt19, Lct, Mdh1, Mrpl12, Muc13, Myl6, Myo15b, Ndufa12, Ndufa7, Ndufs6, Ndufv3, Neat1, Oaz1, Pls1, Ppia, Scp2, Sectm1b, Sepp1, Slc5a1, Smim24, Tceb2, Tmsb10, Tomm7, Tpt1, Ubl5, Uqcr10, Uqcr11, Uqcrb, Uqcrq, Vil1
genesDE <- rownames(DE)[abs(DE$avg_log2FC) > 1 & DE$p_val_adj < 0.05]
# 2010106E10Rik, Abcc2, Adh6a, Agr2, Akr1b8, Aldh1a1, Aldoa, Aldob, Angptl4, Anxa2, Anxa4, Anxa5, Apoa1, Apoa4, Apob, Apoc3, B2m, Calml4, Capg, Carhsp1, Ccl25, Ccl5, Cd7, Cdhr2, Ces1f, Ces2e, Cgref1, Ckb, Clec2h, Crip1, Cyb5r3, Cycs, Cyp27a1, Cyp2b10, Cyp2d26, Cyp3a11, Cyp3a25, Cystm1, Dgat2, Dmbt1, Dstn, Eps8l3, Espn, Ethe1, Exoc3l4, Fabp1, Fam101b, Fam162a, Fkbp5, Flnb, Fth1, Gale, Ggt1, Gjb1, Gna11, Gp1bb, Gpx4, Gsdmd, Gsta4, Gstk1, Gstm1, Gstm3, Guca2a, Guca2b, Gzma, Gzmb, H2-K1, H2-Q1, H2-Q2, H2-T3, Hadh, Hmgcs2, Hsd17b11, Hspa8, Ifi27l2b, Khk, Krt18, Krt19, Krt7, Krt8, Lct, Lgals2, Lgals4, Lmna, Ly6d, Malat1, Maoa, Mat2a, Mgst1, Mgst3, Mmp15, Mt1, Mt2, Mttp, Muc13, Muc3, Myl7, Myo15b, Naprt, Ndrg1, Npc1l1, P4hb, Papss2, Pcyt2, Pepd, Pigr, Plac8, Plin2, Pls1, Plscr1, Prdx1, Qsox1, Rbp2, Reep6, Reg1, Rfk, Ripk3, Rpsa, S100a10, S100a11, S100a6, Scp2, Sepp1, Serpinb1a, Serpinb6a, Sfn, Sgk1, Sh3bgrl3, Sis, Slc26a6, Slc27a4, Slc35c2, Slc51a, Slc5a1, Slc6a19, Slc6a20a, Slc6a6, Slc9a3r1, Sod1, Spink4, Sult1b1, Sult1d1, Sult2b1, Tcn2, Tkfc, Tm4sf20, Tm4sf4, Tm4sf5, Tmem120a, Tmsb10, Treh, Tspan8, Xdh

length(intersect(genesDE, genesDR))/length(genesDE)
# % Overlap between DR and DE
# 0.1176471

enrichmentComparison(genesDE, genesDR)
# Shared DROnly DEONly
#   103    156    110

FC <- DE$avg_log2FC
names(FC) <- rownames(DE)

t.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR], alternative = 'greater')
# Welch Two Sample t-test
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# t = 2.9634, df = 47.113, p-value = 0.00238
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.2642159       Inf
# sample estimates:
#   mean of x   mean of y
# 0.58850575 -0.02055576

wilcox.test(FC[names(FC) %in% genesDR], FC[!names(FC) %in% genesDR])
# Wilcoxon rank sum test with continuity correction
#
# data:  FC[names(FC) %in% genesDR] and FC[!names(FC) %in% genesDR]
# W = 32504, p-value = 0.001409
# alternative hypothesis: true location shift is not equal to 0

plotData <- data.frame(G = names(FC), FC = FC, DR = factor(ifelse(names(FC) %in% genesDR, 'Yes', 'No'), levels = c('Yes', 'No')))
plotData$G[plotData$DR == 'No'] <- NA
plotData$G[abs(plotData$FC) < 1] <- NA
pPos <- position_jitter(seed = 2)
png('HNF4A-HNF4G/Results/drde_Hnf4ag.png', width = 1000, height = 2000, res = 300)
ggplot(plotData, aes(DR, FC, label = G)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  geom_jitter(alpha = .3, pch = 16, position = pPos) +
  xlab('scTenifoldKnk\nDifferentially Regulated') +
  ylab(parse(text = 'log[2]~(Fold-Change)~by~MAST')) +
  labs(title = 'Hnf4ag') +
  theme(plot.title = element_text(face = 2)) +
  geom_abline(slope = 0, intercept = c(-1,1), lty = 2, col = 'red') +
  geom_text_repel(position = pPos, min.segment.length = 0)
dev.off()
