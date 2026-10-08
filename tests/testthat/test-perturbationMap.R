# Counts with three gene modules (same generator as test-heat.R)
simulateModules <- function(nGenes = 60, nCells = 600) {
  set.seed(1)
  modules <- matrix(rgamma(3 * nCells, shape = 2), 3, nCells)
  loadings <- matrix(0, nGenes, 3)
  loadings[cbind(seq_len(nGenes), rep(1:3, length.out = nGenes))] <- runif(nGenes, 0.5, 2)
  X <- matrix(rpois(nGenes * nCells, lambda = 5 * loadings %*% modules), nGenes, nCells)
  rownames(X) <- paste0('g', seq_len(nGenes))
  colnames(X) <- paste0('c', seq_len(nCells))
  X
}

O <- suppressMessages(scTenifoldKnk(simulateModules(), transcriptomeWide = TRUE, qc = FALSE,
                                    nc_nNet = 3, nc_nCells = 300, nCores = 1))
profileOf <- function(g) {
  p <- rank(O$perturbationDistances[g, ]) / ncol(O$perturbationDistances) * O$perturbationDirections[g, ]
  p[g] <- 0
  p
}

test_that("the knockout whose profile generated the signature ranks first", {
  set.seed(2)
  sig <- profileOf('g7') + rnorm(60, sd = 0.05)
  M <- perturbationMap(O, signature = sig, genes = 'g7', plot = FALSE)
  expect_equal(M$gene[1], 'g7')
  expect_equal(M$rank[1], 1)
  expect_true(all(diff(M$similarity) <= 0))
  expect_true(all(M$similarity >= -1 - 1e-12 & M$similarity <= 1 + 1e-12))
  expect_equal(M$x, M$similarity)
  expect_named(M, c('gene', 'similarity', 'rank', 'percentile', 'x', 'y'))
})

test_that("similarity is the cosine between the signed profile and the signature", {
  set.seed(3)
  sig <- setNames(rnorm(60), colnames(O$perturbationDistances))
  M <- perturbationMap(O, signature = sig, plot = FALSE)
  p <- profileOf('g20')
  expect_equal(M$similarity[M$gene == 'g20'], sum(p * sig) / sqrt(sum(p^2) * sum(sig^2)))
  # a reversed signature reverses every similarity
  R <- perturbationMap(O, signature = -sig, plot = FALSE)
  expect_equal(R$similarity[match(M$gene, R$gene)], -M$similarity)
})

test_that("the PCA layout projects the signature and both layouts plot", {
  sig <- profileOf('g7')
  M <- perturbationMap(O, signature = sig, layout = 'pca', plot = FALSE)
  expect_length(attr(M, 'signature'), 2)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(perturbationMap(O, signature = sig, genes = c('g7', 'g8'), layout = 'pca'))
  expect_silent(perturbationMap(O, signature = sig, genes = 'g7'))
})

test_that("invalid inputs are reported", {
  sig <- profileOf('g7')
  expect_error(perturbationMap(O[c('tensorNetworks', 'perturbationDistances')], signature = sig),
               'dr_direction')
  expect_error(perturbationMap(O, signature = unname(sig)), 'named numeric')
  expect_error(perturbationMap(O, signature = c(x = 1, y = 2)), 'fewer than 10')
  expect_warning(perturbationMap(O, signature = sig, genes = 'notAGene', plot = FALSE), 'notAGene')
})
