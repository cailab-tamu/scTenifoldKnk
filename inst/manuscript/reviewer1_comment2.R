library(igraph)

# NKX2-1
load('NKX2-1/Results/GSM3716703.RData')
N <- graph_from_adjacency_matrix(GSM3716703$WT, weighted = TRUE)
D <- degree(N)
NKX21 <- D['Nkx2-1']

# TREM2
load('TREM2/Results/GSE130626.RData')
N <- graph_from_adjacency_matrix(GSE130626$WT, weighted = TRUE)
D <- degree(N)
TREM2 <- D['Trem2']

# HNF4AG
load('HNF4A-HNF4G/Results/GSM3477499.RData')
N <- graph_from_adjacency_matrix(GSM3477499$WT, weighted = TRUE)
D <- degree(N)
HNF4AG <- D[c('Hnf4a','Hnf4g')]

# CFTR
load('CFTR/Results/SRS4245406.RData')
N <- graph_from_adjacency_matrix(SRS4245406$WT, weighted = TRUE)
D <- degree(N)
CFTR <- D['Cftr']

# DMD
load('DMD/Results/GSM4116571.RData')
N <- graph_from_adjacency_matrix(GSM4116571$WT, weighted = TRUE)
D <- degree(N)
DMD <- D['Dmd']

# MECP2r1
load('MECP2/Results/SRS3059998.RData')
N <- graph_from_adjacency_matrix(SRS3059998$WT, weighted = TRUE)
D <- degree(N)
MECP2r1 <- D['Mecp2']

# MECP2r1
load('MECP2/Results/SRS3059999.RData')
N <- graph_from_adjacency_matrix(SRS3059999$WT, weighted = TRUE)
D <- degree(N)
MECP2r2 <- D['Mecp2']

mean(c(NKX21, TREM2, HNF4AG, CFTR, DMD, MECP2r1, MECP2r2))
# Average degree of reported knockouts
# 2475.125
sd(c(NKX21, TREM2, HNF4AG, CFTR, DMD, MECP2r1, MECP2r2))
# Standard deviation of the degree for the reported knockouts
# 2475.125