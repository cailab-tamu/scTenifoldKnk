library(pbapply)
library(ggplot2)
library(ggridges)

fileList <- list.files(pattern = '.RData')

P <- lapply(1:10, function(X){
  # load(fileList[X])
  # DR <- O$diffRegulation
  # write.csv(DR, file = paste0('r',X,'sko_Trem2.csv'))
  DR <- read.csv(file = paste0('r',X,'sko_Trem2.csv'), row.names = 1)
  Z <- DR$Z
  names(Z) <- DR$gene
  return(Z)
})

geneList <- unique(unlist(lapply(P, names)))

P <- lapply(P, function(X){X[geneList]})
P <- as.data.frame(P)
P <- as.matrix(P)
P <- preprocessCore::normalize.quantiles(P)
rownames(P) <- geneList
colnames(P) <- paste0('R', 1:10)
P <- P[complete.cases(P),]

sKO <- cor(P, method = 'sp')
boxplot(sKO[lower.tri(sKO, diag = FALSE)])

P <- pblapply(1:10, function(X){
  # load(fileList[X])
  # WT <- O$tensorNetworks$WT
  # KO <- O$tensorNetworks$WT
  # KO[c('Trem2','Apoe'),] <- 0
  # MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  # DR <- scTenifoldKnk:::dRegulation(MA, c('Trem2', 'Apoe'))
  # write.csv(DR, file = paste0('r',X,'dko_Trem2-Apoe.csv'))
  DR <- read.csv(file = paste0('r',X,'dko_Trem2-Apoe.csv'), row.names = 1)
  Z <- DR$Z
  names(Z) <- DR$gene
  return(Z)
})

geneList <- unique(unlist(lapply(P, names)))

P <- lapply(P, function(X){X[geneList]})
P <- as.data.frame(P)
P <- as.matrix(P)
P <- preprocessCore::normalize.quantiles(P)
rownames(P) <- geneList
colnames(P) <- paste0('R', 1:10)
P <- P[complete.cases(P),]

dKO <- cor(P, method = 'sp')
boxplot(dKO[lower.tri(dKO, diag = FALSE)])


P <- pblapply(1:10, function(X){
  # load(fileList[X])
  # WT <- O$tensorNetworks$WT
  # KO <- O$tensorNetworks$WT
  # KO[c('Trem2', 'Apoe', 'Lpl'),] <- 0
  # MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  # DR <- scTenifoldKnk:::dRegulation(MA, c('Trem2', 'Apoe', 'Lpl'))
  # write.csv(DR, file = paste0('r',X,'tko_Trem2-Apoe-Lpl.csv'))
  DR <- read.csv(file = paste0('r',X,'tko_Trem2-Apoe-Lpl.csv'), row.names = 1)
  Z <- DR$Z
  names(Z) <- DR$gene
  return(Z)
})

geneList <- unique(unlist(lapply(P, names)))

P <- lapply(P, function(X){X[geneList]})
P <- as.data.frame(P)
P <- as.matrix(P)
P <- preprocessCore::normalize.quantiles(P)
rownames(P) <- geneList
colnames(P) <- paste0('R', 1:10)
P <- P[complete.cases(P),]

tKO <- cor(P, method = 'sp')
boxplot(tKO[lower.tri(tKO, diag = FALSE)])

corValues <- data.frame(sKO[lower.tri(sKO, diag = FALSE)], dKO[lower.tri(dKO, diag = FALSE)], tKO[lower.tri(tKO, diag = FALSE)])
colnames(corValues) <- c('Trem2', 'Trem2 + Apoe', 'Trem2 + Apoe + Lpl')
corValues <- reshape2::melt(corValues)
colnames(corValues) <- c('KO', 'rho')
rownames(corValues) <- NULL


png('stabilityResults.png', width = 1500, height = 750, res = 300)
ggplot(corValues, aes(x = rho, y = KO)) +
  geom_density_ridges(jittered_points = TRUE,
                      quantile_lines = TRUE,
                      scale = 0.6, vline_size = 1, vline_color = "red", alpha = 0.5,
  position = position_raincloud(adjust_vlines = TRUE)
) + theme_ridges() + ylab('Knockout') + xlab(parse(text = 'rho'))
dev.off()


corValues <- split(corValues, corValues$KO)
paste0(round(unlist(lapply(corValues, function(X){mean(X[,2])})),2), ' +/- ',
round(unlist(lapply(corValues, function(X){sd(X[,2])})),2))

