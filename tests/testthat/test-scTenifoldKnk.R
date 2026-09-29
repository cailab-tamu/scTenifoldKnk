# Counts with gene modules so that the networks have strong edges
simulateCounts <- function(nGenes = 60, nCells = 600) {
  set.seed(1)
  modules <- matrix(rgamma(3 * nCells, shape = 2), 3, nCells)
  loadings <- matrix(0, nGenes, 3)
  loadings[cbind(seq_len(nGenes), rep(1:3, length.out = nGenes))] <- runif(nGenes, 0.5, 2)
  X <- matrix(rpois(nGenes * nCells, lambda = 5 * loadings %*% modules), nGenes, nCells)
  rownames(X) <- paste0('g', seq_len(nGenes))
  colnames(X) <- paste0('c', seq_len(nCells))
  X
}

runKnk <- function(X, ...) {
  suppressMessages(scTenifoldKnk(X, qc = FALSE, nc_nNet = 3, nc_nCells = 300,
                                 nCores = 1, ...))
}

test_that("single-gene knockout: the gene is left out of the expectation", {
  O <- runKnk(simulateCounts(), gKO = 'g1')
  DR <- O$diffRegulation
  expect_equal(DR$gene[1], 'g1')
  expect_equal(mean(DR$FC[DR$gene != 'g1']), 1)
  expect_equal(DR, dRegulation(O$manifoldAlignment, gKO = 'g1'))
})

test_that("only the WT network is returned", {
  O <- runKnk(simulateCounts(), gKO = c('g1', 'g2'))
  expect_named(O, c('tensorNetworks', 'manifoldAlignment', 'diffRegulation'))
  expect_named(O$tensorNetworks, 'WT')
  # The KO network used in the alignment is WT with the gKO rows set to 0
  KO <- as.matrix(O$tensorNetworks$WT)
  KO[c('g1', 'g2'), ] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(as.matrix(O$tensorNetworks$WT), KO, d = 2, nCores = 1)
  expect_equal(MA, O$manifoldAlignment)
})

test_that("multi-gene knockout: all genes are left out of the expectation", {
  gKO <- c('g1', 'g2', 'g4')
  O <- runKnk(simulateCounts(), gKO = gKO)
  DR <- O$diffRegulation
  expect_equal(mean(DR$FC[!DR$gene %in% gKO]), 1)
  expect_true(all(gKO %in% DR$gene))
  expect_equal(DR, dRegulation(O$manifoldAlignment, gKO = gKO))
})

test_that("transcriptome-wide mode still reports the distances", {
  O <- runKnk(simulateCounts(), gKO = c('g1', 'g2'), transcriptomeWide = TRUE)
  expect_equal(dim(O$perturbationDistances), c(2, 60))
  expect_false(anyNA(O$perturbationDistances))
})
