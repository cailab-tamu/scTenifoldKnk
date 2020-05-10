library(Seurat)
library(Matrix)
load('C:/Users/danie/Filr/Net Folders/VIBS Lab/Cailab/ChapkinLab_CailabReport/CailabReport_Dataset.RData')
metadata <- read.csv('C:/Users/danie/Filr/Net Folders/VIBS Lab/Cailab/ChapkinLab_CailabReport/CailabReport_Table2.csv', stringsAsFactors = FALSE)
Idents(ALL) <- paste0(metadata$CT, '_', metadata$CLASS)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)

sapply(unique(metadata$CT), function(X){
  message(X)
  DE <- FindMarkers(ALL, ident.1 = paste0(X,'_WT'), ident.2 = paste0(X,'_KO'), test.use = 'MAST', logfc.threshold = 0)
  write.csv(DE, file = paste0('DE_',X,'.csv'))
})
