library(UpSetR)
library(pbapply)
library(enrichR)


folder <- 'all'
q <- 0.8

load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'CFTR_',q,'_SRS3161261.RData'))
CFTR_1 <- SRS3161261
rm(SRS3161261)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'CFTR_',q,'_SRS4245406.RData'))
CFTR_2 <- SRS4245406
rm(SRS4245406)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'HNF4AG_',q,'_GSM4116571.RData'))
HNF4AG <- GSM3477499
rm(GSM3477499)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'HNF4ASMAD4_',q,'_GSM4116571.RData'))
HNF4ASMAD4 <- GSM3477499
rm(GSM3477499)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'MECP2_',q,'_SRS3059998.RData'))
MECP2_1 <- WT
rm(WT)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'MECP2_',q,'_SRS3059999.RData'))
MECP2_2 <- WT
rm(WT)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'NKX21_',q,'_GSM3716703.RData'))
NKX21 <- WT
rm(WT)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'TREM2_',q,'_GSE130626.RData'))
TREM2 <- WT
rm(WT)
load(paste0(folder,q,'/',ifelse(folder=='all','','nmr'),'DMD_',q,'_GSM4116571.RData'))
DMD <- GSM4116571
rm(GSM4116571)

gList_CF <- readLines('gSet_CF.txt')
gList_DMD <- readLines('gSet_DMD.txt')
gList_MECP2 <- readLines('gSet_Rett.txt')
gList_HNF4AG <- c('SPDEF', 'ATOH1','MUC2','TFF3')

allDatasets <- list(CFTR_1, CFTR_2, HNF4AG, HNF4ASMAD4, MECP2_1, MECP2_2, NKX21, TREM2, DMD)
nGenes <- pbsapply(seq_len(30)[-1], function(Z){
  dGenes <- lapply(allDatasets, function(X){
    X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment[,seq_len(Z)])
    X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
    if('Hnf4a' %in% X$diffRegulation$gene[1:2]){
      D <- X$diffRegulation[-c(1,2),]$distance
      PV <- c(0,0,pchisq((D^2/mean(D^2)), df = 1, lower.tail = FALSE))
    } else{
      D <- X$diffRegulation[-1,]$distance
      PV <- c(0,pchisq((D^2/mean(D^2)), df = 1, lower.tail = FALSE))
    }
    PA <- p.adjust(PV, method = 'fdr')
    dGenes <- X$diffRegulation$gene[PA < 0.05]
  })
  M1 <- mean(toupper(dGenes[[1]]) %in% gList_CF)
  M2 <- mean(toupper(dGenes[[2]]) %in% gList_CF)
  M3 <- mean(toupper(dGenes[[3]]) %in% gList_HNF4AG)
  M5 <- mean(toupper(dGenes[[5]]) %in% gList_MECP2)
  M6 <- mean(toupper(dGenes[[6]]) %in% gList_MECP2)
  M9 <- mean(toupper(dGenes[[9]]) %in% gList_DMD)
  M10 <- mean(rowSums(fromList(dGenes[1:2])) == 2)
  M11 <- mean(rowSums(fromList(dGenes[3:4])) == 2)
  M12 <- mean(rowSums(fromList(dGenes[5:6])) == 2)
  c(M1,M2,M3,M5,M6,M9,M10,M11,M12)
})
Z <- which.max(rowMeans(apply(nGenes,1,rank))) + 1


dbList <- c("BioPlanet_2019", "KEGG_2019_Mouse", "Reactome_2016","GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
dGenes <- lapply(allDatasets, function(X){
  X$diffRegulation <- scTenifoldNet::dRegulation(X$manifoldAlignment[,seq_len(2)])
  X$diffRegulation <- X$diffRegulation[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^mt-', X$diffRegulation$gene, ignore.case = TRUE),]
  if('Hnf4a' %in% X$diffRegulation$gene[1:2]){
    D <- X$diffRegulation[-c(1,2),]$distance
    PV <- c(0,0,pchisq((D^2/mean(D^2)), df = 1, lower.tail = FALSE))
  } else{
    D <- X$diffRegulation[-1,]$distance
    PV <- c(0,pchisq((D^2/mean(D^2)), df = 1, lower.tail = FALSE))
  }
  PA <- p.adjust(PV, method = 'fdr')
  dGenes <- X$diffRegulation$gene[PA < 0.05]
  return(dGenes)
  # E <- enrichr(dGenes, dbList)
  # E <- do.call(rbind.data.frame,E)
  # E <- E[E$Adjusted.P.value < 0.05,]
  # E <- E[order(E$P.value),c(1)]
  # return(E)
})
names(dGenes) <- c("CFTR_1", "CFTR_2", "HNF4AG", "HNF4ASMAD4", "MECP2_1", "MECP2_2", "NKX21", "TREM2", 'DMD')
lengths(dGenes)

upset(fromList(dGenes), nsets = length(dGenes))

library(fgsea)

