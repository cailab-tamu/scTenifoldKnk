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
| 3 | `makeNetworks` | Constructs gene regulatory networks from subsampled cells using principal component regression (`pcNet`) for both WT and KO conditions |
| 4 | `tensorDecomposition` | CANDECOMP/PARAFAC (CP) tensor decomposition for network denoising |
| 5 | `manifoldAlignment` | Non-linear manifold alignment of the WT and KO denoised networks |
| 6 | `dRegulation` | Differential regulation testing via Box-Cox transformation and chi-square statistics |

Individual functions are exported and fully documented, allowing users to run or modify any step independently.

## Input

The required input is a **raw counts matrix** with genes as rows and cells (barcodes) as columns. Data should be *unnormalized* when `qc = TRUE` (the default). The modular design allows users to substitute custom preprocessing at any step.

## Output

`scTenifoldKnk()` returns a list with three elements:

- **`tensorNetworks`** — Weight-averaged denoised gene regulatory networks after CP tensor decomposition, containing:
  - `WT`: The network for the wild-type condition (sparse matrix of class `dgCMatrix`).
  - `KO`: The network for the knocked-out condition (sparse matrix of class `dgCMatrix`).
- **`manifoldAlignment`** — A data frame of low-dimensional features from the non-linear manifold alignment, with 2 × *n* genes rows and *d* columns (default *d* = 2).
- **`diffRegulation`** — A data frame with six columns:
  - `gene`: Gene identifier.
  - `distance`: Euclidean distance between the gene's coordinates in the two conditions.
  - `Z`: Z-score after Box-Cox power transformation.
  - `FC`: Fold change with respect to the expectation.
  - `p.value`: P-value from the chi-square distribution with one degree of freedom.
  - `p.adj`: Adjusted p-value (Benjamini & Hochberg FDR correction).

## Running Time

Running time is largely determined by the network construction step and scales with the number of cells and genes. Representative benchmarks:

| Cells | Genes | Time |
|------:|------:|-----:|
| 300 | 1,000 | 3.45 min |
| 1,000 | 1,000 | 4.25 min |
| 1,000 | 5,000 | 2 h 51.6 min |
| 2,500 | 5,000 | 2 h 55.3 min |
| 5,000 | 5,000 | 3 h 8.9 min |
| 5,000 | 7,500 | 3 h 9.5 min |
| 7,500 | 5,000 | 10 h 15.5 min |
| 7,500 | 7,500 | 10 h 16.1 min |

## Example

### Simulating a dataset

We create a sparse count matrix of 2,000 cells and 100 genes drawn from a negative binomial distribution (~67 % zeros). The last ten genes are prefixed with `mt-` to simulate mitochondrial genes.

```r
library(scTenifoldKnk)

nCells <- 2000
nGenes <- 100
set.seed(1)
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

See also: [plotKO() — Frequently Asked Questions](plotKO_FAQ.md)

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
