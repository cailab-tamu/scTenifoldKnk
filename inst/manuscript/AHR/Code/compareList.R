compareList <- function(gList1, gList2, B = 100, label = 'L'){
  require(pbapply)
  require(ggplot2)
  overlapR <- pbsapply(seq_len(B), function(Z){
    nRanks <- min(c(length(gList1), length(gList2)))
    rList1 <- sample(gList1)
    rList2 <- sample(gList2)
    rOverlap <- sapply(seq_len(nRanks), function(X){
      length(intersect(rList1[seq_len(X)], rList2[seq_len(X)]))
    })
    return(rOverlap)
  })
  
  nRanks <- min(c(length(gList1), length(gList2)))
  rOverlap <- sapply(seq_len(nRanks), function(X){
    length(intersect(gList1[seq_len(X)], gList2[seq_len(X)]))
  })
  
  listComparison <- data.frame(
    rank = seq_len(nRanks),
    real = rOverlap,
    expLB = apply(overlapR,1,mean) - 5* apply(overlapR,1,sd),
    exp = apply(overlapR,1,mean),
    expUB = apply(overlapR,1,mean) + 5* apply(overlapR,1,sd)
  )
  #png('../Results/cl_MECP2.png', width = 1800, height = 1000, res = 300, pointsize = 20)
  ggplot(data = listComparison, mapping = aes(x = rank, y = real)) +
    geom_line(mapping = aes(colour='WT\nKO')) +
    geom_ribbon(mapping = aes(ymin=expLB, ymax=expUB), alpha=0.3, fill = "#E69F00") +
    geom_line(mapping = aes(x = rank, y = exp, colour='Expected by\nrandom')) +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    theme_bw() + xlab('Rank') + ylab('Size of Overlap') + ylim(c(0,7776)) + xlim(c(0,7776)) +
    scale_fill_manual(name="Overlap:",aesthetics = 'colour', values=c("#E69F00","#999999")) +
    theme(legend.position="top")
  #dev.off()  
}
