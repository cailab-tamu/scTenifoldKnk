# scTenifoldKnk

[![CRAN](https://www.r-pkg.org/badges/version/scTenifoldKnk)](https://CRAN.R-project.org/package=scTenifoldKnk)
[![License: GPL (>=2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D2%29-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

**scTenifoldKnk** is an R package for performing virtual knockout experiments on single-cell gene regulatory networks (scGRNs). It uses single-cell RNA-seq (scRNA-seq) data from wild-type (WT) control samples to construct a scGRN, then simulates a gene knockout by zeroing the target gene's outdegree edges in the adjacency matrix. The resulting knocked-out scGRN is compared with the WT scGRN to identify differentially regulated genes, or virtual-knockout perturbed genes, which reveal the functional impact of the knocked-out gene in the analyzed cell population.

Implementations in other languages are also available:

- **Python**: [scTenifoldpy](https://github.com/qwerty239qwe/scTenifoldpy)
- **MATLAB**: [scGEAToolbox](https://github.com/jamesjcai/scGEAToolbox)

## Installation

**scTenifoldKnk** is available on CRAN:

```r
install.packages("scTenifoldKnk")
```

To install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("cailab-tamu/scTenifoldKnk")
```

## Pipeline Overview

The `scTenifoldKnk()` function orchestrates a virtual knockout pipeline built on top of the **scTenifoldNet** framework. Each step reports progress to the console via the [cli](https://cli.r-lib.org/) package.

| Step | Function | Description |
|:----:|:---------|:------------|
| 1 | `scQC` | Quality control — filters cells by library size, outlier detection, minimum gene expression fraction, and mitochondrial read ratio |
| 2 | `cpmNormalization` | Counts-per-million (CPM) normalization |
| 3 | `makeNetworks` | Constructs gene regulatory networks from subsampled cells using principal component regression (`pcNet`). All the per-gene regressions come from one eigendecomposition, which gives the same networks as fitting each gene separately |
| 4 | `tensorDecomposition` | CANDECOMP/PARAFAC (CP) tensor decomposition for network denoising |
| 5 | `strictDirection` | Enforces directionality of the reconstructed adjacency matrix |
| 6 | `manifoldAlignment` | Non-linear manifold alignment of the WT and KO denoised networks |
| 7 | `dRegulation` | Differential regulation testing via Box-Cox transformation and chi-square statistics |

Individual functions are exported and fully documented, allowing users to run or modify any step independently.

## Input

The required input is a **raw counts matrix** with genes as rows and cells (barcodes) as columns. Data should be *unnormalized* when `qc = TRUE` (the default). The modular design allows users to substitute custom preprocessing at any step.

## Operating Modes

`scTenifoldKnk()` supports two modes, selected with the `transcriptomeWide` argument:

- **Knockout** (`transcriptomeWide = FALSE`, default) — Knocks out one target gene, or several genes together (a character vector in `gKO`, e.g. `c("Hnf4a", "Hnf4g")`), and returns the WT/KO networks, the manifold alignment, and the differential regulation table. If none of the target genes has outgoing edges in the WT network, the knockout leaves the network unchanged and a warning says the results are only numerical noise. Target genes without outgoing edges are also reported when only some of them lack edges.
- **Transcriptome-wide perturbation** (`transcriptomeWide = TRUE`) — Builds the WT network **once**, then knocks out every gene in the network (or a user-supplied subset passed through `gKO`), running the manifold alignment and distance calculation for each. It returns a matrix of perturbation distances instead of a single differential regulation table. Because it performs one alignment per perturbed gene, running time scales with the number of genes.

## Reproducibility

Network construction subsamples cells at random, so results depend on the random seed. The `seed` argument (default `1`) is set before each random step. The same input and parameters therefore give the same result on every run, whatever the caller's RNG state, and that state is restored when the function returns.

- Use different values (`seed = 1`, `seed = 2`, ...) to see how much results vary between runs.
- Use `seed = NULL` to let a `set.seed()` call made before `scTenifoldKnk()` control the result.

## Function Reference

All functions below are exported and individually documented (`?functionName`).

### `scTenifoldKnk()`

Main entry point running the full virtual knockout pipeline.

| Argument | Default | Description |
|:---------|:--------|:------------|
| `countMatrix` | — | Raw counts matrix, genes (symbols) as rows, cells as columns. |
| `gKO` | `NULL` | Knockout mode: gene symbol to knock out, or a character vector of genes to knock out together. Transcriptome-wide mode: optional character vector of genes to perturb, each one separately; `NULL` perturbs every gene in the WT network. |
| `transcriptomeWide` | `FALSE` | If `TRUE`, perturb each target gene in turn and return a distance matrix. |
| `qc` | `TRUE` | Apply quality control (`scQC`) to the input matrix. |
| `qc_minLibSize` | `1000` | Minimum library size for a cell to be retained. |
| `qc_removeOutlierCells` | `TRUE` | Remove cells whose library size is an outlier. |
| `qc_minPCT` | `0.05` | Minimum fraction of cells in which a gene must be expressed. |
| `qc_maxMTratio` | `0.1` | Maximum mitochondrial read ratio per cell (genes matching `^MT-`, case-insensitive). |
| `nc_lambda` | `0` | Directionality weighting applied to the weaker edge between two genes. |
| `nc_nNet` | `10` | Number of PC-regression networks to generate. |
| `nc_nCells` | `500` | Number of cells subsampled per network. |
| `nc_nComp` | `3` | Number of principal components used to build networks. |
| `nc_scaleScores` | `TRUE` | Normalize weights so the maximum absolute value is 1. |
| `nc_symmetric` | `FALSE` | Return a symmetric weights matrix. |
| `nc_q` | `0.9` | Cut-off quantile of top relationships to keep. |
| `nc_priorNetwork` | `NULL` | Optional prior network (`data.frame` with `regulators`/`targets`). |
| `td_K` | `3` | Number of rank-one tensors for CP tensor decomposition. |
| `td_maxIter` | `1000` | Maximum tensor-decomposition iterations. |
| `td_maxError` | `1e-05` | Relative Frobenius-norm error tolerance. |
| `td_nDecimal` | `3` | Number of decimal places retained. |
| `ma_nDim` | `2` | Number of manifold-alignment dimensions. |
| `dr_empiricalNull` | `FALSE` | Assign differential regulation p-values using Efron's empirical null (via `locfdr`) instead of the theoretical chi-square null. |
| `nCores` | `parallel::detectCores()` | Number of cores used by the manifold alignment. |
| `seed` | `1` | Seed set before each random stage, so results are reproducible and do not depend on the caller's RNG state (which is restored on exit). Use different values to assess run-to-run variability, or `NULL` to use the caller's RNG state (e.g. a previous `set.seed()`). |

### `scQC()`

Standalone single-cell quality control. Arguments: `X` (raw counts matrix), `minLibSize`, `removeOutlierCells`, `minPCT`, `maxMTratio`, and an optional `label` for progress messages. Returns a `dgCMatrix` containing the cells and genes that pass the filters.

### `dRegulation()`

Differential regulation testing from a manifold alignment. Arguments: `manifoldOutput` (the labeled `manifoldAlignment` matrix, `X_` genes followed by `Y_` genes in the same order) and `empiricalNull` (if `TRUE`, estimate the null distribution of the Z-scores with Efron's empirical null via the `locfdr` package instead of the theoretical chi-square null). Returns the six-column differential regulation table described under [Output](#output).

### `plotKO()`

Plots the KO-centered subnetwork from a `scTenifoldKnk()` result.

| Argument | Default | Description |
|:---------|:--------|:------------|
| `X` | — | Output list from `scTenifoldKnk()`. |
| `gKO` | — | Gene symbol(s) of the simulated knockout, as passed to `scTenifoldKnk()`. |
| `q` | `0.99` | Edge-weight quantile used to threshold weak edges. |
| `annotate` | `TRUE` | Query enrichment databases (`enrichR`) and overlay category pies on nodes. |
| `nCategories` | `20` | Maximum number of enrichment categories shown in the legend. |
| `fdrThreshold` | `0.05` | Adjusted p-value cutoff for reporting enriched terms. |

See also: [plotKO() — Frequently Asked Questions](plotKO_FAQ.md)

## Output

### Single knockout mode (`transcriptomeWide = FALSE`)

`scTenifoldKnk()` returns a list with three elements:

- **`tensorNetworks`** — Weight-averaged denoised gene regulatory networks after CP tensor decomposition, containing:
  - `WT`: The network for the wild-type condition (a `Matrix` object).
  - `KO`: The network for the knocked-out condition (a `Matrix` object).
- **`manifoldAlignment`** — A data frame of low-dimensional features from the non-linear manifold alignment, with 2 × *n* genes rows and *d* columns (default *d* = 2).
- **`diffRegulation`** — A data frame with six columns:
  - `gene`: Gene identifier.
  - `distance`: Euclidean distance between the gene's coordinates in the two conditions.
  - `Z`: Z-score after Box-Cox power transformation.
  - `FC`: Fold change with respect to the expectation.
  - `p.value`: P-value from the chi-square distribution with one degree of freedom, or from Efron's empirical null when `dr_empiricalNull = TRUE`.
  - `p.adj`: Adjusted p-value (Benjamini & Hochberg FDR correction).

### Transcriptome-wide mode (`transcriptomeWide = TRUE`)

`scTenifoldKnk()` returns a list with two elements:

- **`tensorNetworks`** — A list with the WT weight-averaged denoised gene regulatory network (`WT`).
- **`perturbationDistances`** — A numeric matrix of manifold-alignment distances. Rows are the perturbed genes, columns are all genes in the WT network, and each entry is the distance of a gene under the corresponding knockout.

## Running Time

Running time grows mainly with the number of genes. The number of cells matters little, because each network is built from a fixed-size subsample of cells (`nc_nCells`). Single knockout benchmarks with the default parameters (10 networks of 500 cells) on simulated counts, measured on an Apple M2 Pro (16 GB RAM) with R 4.5 and its reference BLAS. Memory is peak resident memory.

| Cells | Genes | Time | Memory |
|------:|------:|-----:|-------:|
| 300 | 1,000 | 17 s | 1.6 GB |
| 1,000 | 1,000 | 17 s | 1.7 GB |
| 1,000 | 5,000 | 4.0 min | 7.0 GB |
| 2,500 | 5,000 | 3.5 min | 6.4 GB |
| 5,000 | 5,000 | 3.7 min | 5.0 GB |

## Example

### Simulating a dataset

We create a sparse count matrix of 2,000 cells and 100 genes drawn from a negative binomial distribution (~67 % zeros). The last ten genes are prefixed with `mt-` to simulate mitochondrial genes.

```r
library(scTenifoldKnk)

nCells <- 2000
nGenes <- 100
set.seed(1) # seeds the simulated counts; the pipeline seed is set with `seed`
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- round(X)
X <- matrix(X, ncol = nCells)
rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
```

### Running the virtual knockout

```r
output <- scTenifoldKnk(
  countMatrix   = X,
  gKO           = "ng10",
  nc_nNet       = 10,
  nc_nCells     = 500,
  td_K          = 3,
  qc_minLibSize = 30
)
```

### Exploring the output

```r
# Structure of the output
str(output)

# Accessing the WT and KO gene regulatory networks
dim(output$tensorNetworks$WT)
dim(output$tensorNetworks$KO)

# Accessing the manifold alignment result
head(output$manifoldAlignment)

# Differential regulation results — top perturbed genes
head(output$diffRegulation, n = 10)

# Plotting the KO-centered subnetwork
plotKO(output, gKO = "ng10")
```

### Multi-gene knockout

Pass several genes to `gKO` to knock them out together in one simulated experiment:

```r
dkoOutput <- scTenifoldKnk(
  countMatrix   = X,
  gKO           = c("ng10", "ng20"),
  nc_nNet       = 10,
  nc_nCells     = 500,
  td_K          = 3,
  qc_minLibSize = 30
)

head(dkoOutput$diffRegulation, n = 10)
plotKO(dkoOutput, gKO = c("ng10", "ng20"))
```

### Transcriptome-wide perturbation

Knock out every gene in the WT network and collect the manifold-alignment distances for each perturbation:

```r
twOutput <- scTenifoldKnk(
  countMatrix       = X,
  transcriptomeWide = TRUE,
  nc_nNet           = 10,
  nc_nCells         = 500,
  td_K              = 3,
  qc_minLibSize     = 30
)

# Distance matrix: perturbed genes (rows) by all genes (columns)
dim(twOutput$perturbationDistances)
twOutput$perturbationDistances[1:5, 1:5]
```

To restrict the perturbation to a subset of genes, pass them through `gKO`. Each gene is knocked out separately, unlike the multi-gene knockout above:

```r
subset <- scTenifoldKnk(
  countMatrix       = X,
  gKO               = c("ng10", "ng20"),
  transcriptomeWide = TRUE,
  qc_minLibSize     = 30
)
```

## Citation

Osorio, D., Zhong, Y., Li, G., Xu, Q., Yang, Y., Tian, Y., Chapkin, R., Huang, J. Z., & Cai, J. J. (2022). scTenifoldKnk: An Efficient Virtual Knockout Tool for Gene Function Predictions via Single-Cell Gene Regulatory Network Perturbation. *Patterns*, **3**(3), 100434. [doi:10.1016/j.patter.2022.100434](https://doi.org/10.1016/j.patter.2022.100434)

BibTeX:

```bibtex
@Article{osorio2022sctenifoldknk,
  title   = {scTenifoldKnk: An Efficient Virtual Knockout Tool for Gene Function
             Predictions via Single-Cell Gene Regulatory Network Perturbation},
  author  = {Daniel Osorio and Yan Zhong and Guanxun Li and Qian Xu and
             Yongjian Yang and Yanan Tian and Robert Chapkin and
             Jianhua Z. Huang and James J. Cai},
  journal = {Patterns},
  year    = {2022},
  volume  = {3},
  number  = {3},
  pages   = {100434},
  issn    = {2666-3899},
  doi     = {10.1016/j.patter.2022.100434},
}
```

---

&copy; The Texas A&M University System. All rights reserved.
