library(igraph)
library(ggplot2)
library(ggrepel)
library(statsExpressions)

# NKX2-1
load('NKX2-1/Results/GSM3716703.RData')
N <- graph_from_adjacency_matrix(GSM3716703$WT, weighted = TRUE)
D <- degree(N)
dNKX21 <- D['Nkx2-1']
gNKX21 <- length(GSM3716703$diffRegulation$gene[GSM3716703$diffRegulation$p.adj < 0.05])


# TREM2
load('TREM2/Results/GSE130626.RData')
N <- graph_from_adjacency_matrix(GSE130626$WT, weighted = TRUE)
D <- degree(N)
dTREM2 <- D['Trem2']
gTREM2 <- length(GSE130626$diffRegulation$gene[GSE130626$diffRegulation$p.adj < 0.05])

# HNF4AG
load('HNF4A-HNF4G/Results/GSM3477499.RData')
N <- graph_from_adjacency_matrix(GSM3477499$WT, weighted = TRUE)
D <- degree(N)
dHNF4AG <- D[c('Hnf4a','Hnf4g')]
gHNF4AG <- rep(length(GSM3477499$diffRegulation$gene[GSM3477499$diffRegulation$p.adj < 0.05]), 2)

# CFTR
load('CFTR/Results/SRS4245406.RData')
N <- graph_from_adjacency_matrix(SRS4245406$WT, weighted = TRUE)
D <- degree(N)
dCFTR <- D['Cftr']
gCFTR <- length(SRS4245406$diffRegulation$gene[SRS4245406$diffRegulation$p.adj < 0.05])


# DMD
load('DMD/Results/GSM4116571.RData')
N <- graph_from_adjacency_matrix(GSM4116571$WT, weighted = TRUE)
D <- degree(N)
dDMD <- D['Dmd']
gDMD <- length(GSM4116571$diffRegulation$gene[GSM4116571$diffRegulation$p.adj < 0.05])

# MECP2r1
load('MECP2/Results/SRS3059998.RData')
N <- graph_from_adjacency_matrix(SRS3059998$WT, weighted = TRUE)
D <- degree(N)
dMECP2r1 <- D['Mecp2']
gMECP2r1 <- length(SRS3059998$diffRegulation$gene[SRS3059998$diffRegulation$p.adj < 0.05])

# MECP2r1
load('MECP2/Results/SRS3059999.RData')
N <- graph_from_adjacency_matrix(SRS3059999$WT, weighted = TRUE)
D <- degree(N)
dMECP2r2 <- D['Mecp2']
gMECP2r2 <- length(SRS3059999$diffRegulation$gene[SRS3059999$diffRegulation$p.adj < 0.05])


DF <- data.frame(degree = c(dNKX21, dTREM2, dHNF4AG, dCFTR, dDMD, dMECP2r1, dMECP2r2),
                 genes = c(gNKX21, gTREM2, gHNF4AG, gCFTR, gDMD, gMECP2r1, gMECP2r2),
                 id = c('NKX21', 'TREM2', 'HNF4A', 'HNF4G', 'CFTR', 'DMD', 'MECP2r1', 'MECP2r2'))

png('reviewer1_minorComment5.png', width = 1500, height = 1500, res = 300)
ggplot(DF, aes(log(degree), log(genes), label = id)) +
  geom_point() +
  xlab(log~(Network~Degree)) +
  ylab(log~(italic(n)~Significant~Genes~(FDR<0.05))) +
  theme_bw() +
  geom_smooth(method = 'lm', color = 'red', lty = 2, alpha = 0.1) +
  geom_text_repel() +
  labs(subtitle = corr_test(log(DF[,1:2]), x = degree, y = genes, type = 'p')$expression[[1]]) +
  theme(plot.subtitle = element_text(size = 10))
dev.off()

mean(c(NKX21, TREM2, HNF4AG, CFTR, DMD, MECP2r1, MECP2r2))
# Average degree of reported knockouts
# 2475.125
sd(c(NKX21, TREM2, HNF4AG, CFTR, DMD, MECP2r1, MECP2r2))
# Standard deviation of the degree for the reported knockouts
# 2823.621
