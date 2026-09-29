# Manifold alignment output with small distances for most genes and large
# ones for the genes moved in the second condition
simulateManifold <- function(moved = character(0), nGenes = 200, d = 2) {
  set.seed(1)
  geneList <- paste0('g', seq_len(nGenes))
  X <- matrix(rnorm(nGenes * d), nGenes, d)
  Y <- X + matrix(rnorm(nGenes * d, sd = 1e-3), nGenes, d)
  Y[geneList %in% moved, ] <- Y[geneList %in% moved, ] + 1
  M <- rbind(X, Y)
  rownames(M) <- c(paste0('X_', geneList), paste0('Y_', geneList))
  M
}

test_that("the knocked-out gene is left out of the expectation", {
  M <- simulateManifold('g1')
  DR <- dRegulation(M, gKO = 'g1')
  d <- setNames(DR$distance, DR$gene)
  expect_equal(DR$FC, DR$distance^2 / mean(d[names(d) != 'g1']^2))
  expect_equal(mean(DR$FC[DR$gene != 'g1']), 1)
  expect_equal(DR$gene[1], 'g1')
  expect_equal(nrow(DR), 200)
})

test_that("all knocked-out genes are left out of the expectation", {
  gKO <- c('g1', 'g50', 'g100')
  M <- simulateManifold(gKO)
  DR <- dRegulation(M, gKO = gKO)
  d <- setNames(DR$distance, DR$gene)
  expect_equal(DR$FC, DR$distance^2 / mean(d[!names(d) %in% gKO]^2))
  expect_setequal(DR$gene[1:3], gKO)
  expect_true(all(DR$p.adj[DR$gene %in% gKO] < 0.05))
})

test_that("including the knocked-out gene in the expectation hides the other genes", {
  M <- simulateManifold(c('g1', 'g2'))
  # g2 moves a tenth of what g1 moves
  n <- nrow(M) / 2
  M[n + 2, ] <- M[2, ] + (M[n + 2, ] - M[2, ]) / 10
  withKO <- dRegulation(M)
  withoutKO <- dRegulation(M, gKO = 'g1')
  expect_false(withKO$p.adj[withKO$gene == 'g2'] < 0.05)
  expect_true(withoutKO$p.adj[withoutKO$gene == 'g2'] < 0.05)
  # Distances and Z do not depend on gKO
  expect_equal(withKO[order(withKO$gene), c('distance', 'Z')],
               withoutKO[order(withoutKO$gene), c('distance', 'Z')],
               ignore_attr = TRUE)
})

test_that("without gKO the expectation uses all genes", {
  M <- simulateManifold('g1')
  DR <- dRegulation(M)
  expect_equal(DR$FC, DR$distance^2 / mean(DR$distance^2))
})

test_that("invalid gKO values are rejected", {
  M <- simulateManifold('g1', nGenes = 5)
  expect_error(dRegulation(M, gKO = 'absent'), 'not present in the manifold alignment')
  expect_error(dRegulation(M, gKO = TRUE), "must be a character vector")
  expect_error(dRegulation(M, gKO = paste0('g', 1:5)), 'At least one gene')
})
