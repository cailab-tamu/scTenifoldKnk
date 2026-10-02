# Counts with three gene modules (same generator as test-scTenifoldKnk.R)
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
moduleOf <- function(g) ((as.integer(sub('g', '', g)) - 1) %% 3) + 1

runKnk <- function(X, ...) {
  suppressMessages(scTenifoldKnk(X, qc = FALSE, nc_nNet = 3, nc_nCells = 300,
                                 nCores = 1, ...))
}

test_that("heatKernel is the identity at t = 0, symmetric, and matches the eigendecomposition", {
  set.seed(1)
  A <- matrix(rnorm(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  expect_equal(unname(heatKernel(A, t = 0)), diag(6))
  H <- heatKernel(A, t = 3)
  expect_true(isSymmetric(H))
  expect_equal(dimnames(H), dimnames(A))
  S <- (A + t(A)) / 2
  E <- eigen(S, symmetric = TRUE)
  ref <- E$vectors %*% diag(exp(3 * (E$values / max(E$values) - 1))) %*% t(E$vectors)
  expect_equal(unname(H), ref)
  expect_error(heatKernel(A, t = -1))
  expect_error(heatKernel(A, symmetric = FALSE))
})

test_that("knockoutDirection: genes of the knocked-out gene's module go down", {
  X <- simulateModules()
  D <- knockoutDirection(X, gKO = 'g1')
  expect_equal(dim(D), c(1, nrow(X)))
  same <- setdiff(rownames(X)[moduleOf(rownames(X)) == 1], 'g1')
  expect_true(all(D[1, same] < 0))
  # a multi-gene knockout is one row named by the joined genes
  expect_equal(rownames(knockoutDirection(X, gKO = list(c('g1', 'g4')))), 'g1+g4')
})

test_that("default diffRegulation gains direction columns; existing columns are unchanged", {
  X <- simulateModules()
  O <- runKnk(X, gKO = 'g1')
  DR <- O$diffRegulation
  expect_equal(colnames(DR), c('gene', 'distance', 'Z', 'FC', 'p.value', 'p.adj',
                               'direction', 'directionScore'))
  expect_equal(DR[, 1:6], dRegulation(O$manifoldAlignment, gKO = 'g1'))
  expect_equal(DR$direction[DR$gene == 'g1'], 'down')
  expect_true(all(DR$direction %in% c('up', 'down')))
  # dr_direction = FALSE restores the previous output
  O0 <- runKnk(X, gKO = 'g1', dr_direction = FALSE)
  expect_equal(colnames(O0$diffRegulation), c('gene', 'distance', 'Z', 'FC', 'p.value', 'p.adj'))
})

test_that("heat manifold alignment: single knockout ranks the knocked-out gene first", {
  X <- simulateModules()
  O <- runKnk(X, gKO = 'g1', ma_method = 'heat')
  expect_named(O, c('tensorNetworks', 'diffRegulation'))
  expect_equal(O$diffRegulation$gene[1], 'g1')
  expect_equal(mean(O$diffRegulation$FC[O$diffRegulation$gene != 'g1']), 1)
  # same ranking direction as the manifold alignment
  M <- runKnk(X, gKO = 'g1')
  dm <- setNames(M$diffRegulation$distance, M$diffRegulation$gene)
  dh <- setNames(O$diffRegulation$distance, O$diffRegulation$gene)
  g <- setdiff(names(dm), 'g1')
  expect_gt(cor(dm[g], dh[g], method = 'spearman'), 0.3)
})

test_that("transcriptome-wide mode uses the heat kernel by default and reports directions", {
  X <- simulateModules()
  O <- runKnk(X, gKO = c('g1', 'g2'), transcriptomeWide = TRUE)
  expect_named(O, c('tensorNetworks', 'perturbationDistances', 'perturbationDirections'))
  expect_equal(dim(O$perturbationDistances), c(2, 60))
  expect_equal(dim(O$perturbationDirections), c(2, 60))
  expect_true(all(O$perturbationDirections %in% c(-1, 0, 1)))
  # equivalent to calling hkManifoldAlignment on the returned WT network
  expect_equal(O$perturbationDistances,
               hkManifoldAlignment(O$tensorNetworks$WT, X, gKO = c('g1', 'g2')))
  # the manifold alignment remains available
  M <- runKnk(X, gKO = c('g1', 'g2'), transcriptomeWide = TRUE, ma_method = 'manifold')
  expect_false(anyNA(M$perturbationDistances))
})
