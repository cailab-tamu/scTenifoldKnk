# MA <- O$manifoldAlignment
# MA <- MA[!grepl('_MT-|_RPL|_RPS',rownames(MA), ignore.case = TRUE),]
# DR <- scTenifoldNet::dRegulation(MA)
# DR$FC <- (DR$distance^2)/mean(DR$distance[-1]^2)
# DR$p.value <- pchisq(DR$FC, df = 1, lower.tail = FALSE)
# DR$p.adj <- p.adjust(DR$p.value, method = 'fdr')
# writeLines(DR$gene[DR$p.adj < 0.05])
