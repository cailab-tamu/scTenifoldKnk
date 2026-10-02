# scTenifoldKnk

[![R-CMD-check](https://github.com/cailab-tamu/scTenifoldKnk/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cailab-tamu/scTenifoldKnk/actions/workflows/R-CMD-check.yaml)
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
| 6 | `manifoldAlignment` or `hkManifoldAlignment` | Comparison of the WT and KO denoised networks: non-linear manifold alignment (one alignment per knockout), or heat manifold alignment (one heat kernel of the WT network, read for every knockout) |
| 7 | `dRegulation` | Differential regulation testing via Box-Cox transformation and chi-square statistics, with the predicted direction (up/down) of each gene from `knockoutDirection` |

Individual functions are exported and fully documented, allowing users to run or modify any step independently.

## Input

The required input is a **raw counts matrix** with genes as rows and cells (barcodes) as columns, as a `matrix` or a sparse `dgCMatrix`. A data.frame is not accepted: convert it first with `as.matrix()`. Data should be *unnormalized* when `qc = TRUE` (the default). The modular design allows users to substitute custom preprocessing at any step.

## Quality Control: Best Practices

These recommendations come from benchmarking virtual knockouts against bulk knockdown/knockout RNA-seq profiles of the same cell lines (about 1,200 knockout experiments in 17 cell lines, with WT data from three single-cell atlases).

- **Remove ribosomal and mitochondrial genes before building the networks.** Ribosomal protein genes (including their pseudogenes) and mitochondrial genes form dense, highly co-expressed modules that otherwise dominate the networks and crowd the top of the differential regulation table with ribosomal genes, whatever gene is knocked out. Removing them was the largest single improvement found:

  ```r
  riboMito <- grepl('^(RP[LS][0-9]|RPLP[0-9]|RPSA|MRP[LS][0-9]|MT-|MTND[0-9]|MTCO[0-9]|MTATP[0-9]|MTCYB|MTRNR2L)',
                    toupper(rownames(X)))
  X <- X[!riboMito, ]
  ```

  Apply the mitochondrial read filter (`qc_maxMTratio`) before removing these genes, for example by running `scQC()` first.
- **Use about 3,000 highly variable genes** (for example selected with Seurat's `vst` method), always keeping the genes to knock out. Results are stable between 3,000 and 5,000 genes, while running time and memory grow with the square of the number of genes (see [Running Time](#running-time)).
- **Adapt the mitochondrial filter to the data.** The default `qc_maxMTratio = 0.1` suits most data sets, but some platforms (for example 10x v3 libraries of cultured cell lines) have a higher baseline mitochondrial fraction; a per-sample cut-off such as the median plus three median absolute deviations keeps the healthy cells. Mitochondrial genes are detected by gene symbol (`^MT-`), so use gene symbols as row names.
- **Average over seeds when ranking genes for a single knockout matters.** Rankings of the same knockout agree at a Spearman correlation of about 0.9 between seeds; averaging two or three seeds (`seed = 1, 2, 3`) smooths that variation.

## Operating Modes

`scTenifoldKnk()` supports two modes, selected with the `transcriptomeWide` argument:

- **Knockout** (`transcriptomeWide = FALSE`, default) — Knocks out one target gene, or several genes together (a character vector in `gKO`, e.g. `c("Hnf4a", "Hnf4g")`), and returns the WT/KO networks, the manifold alignment, and the differential regulation table. If none of the target genes has outgoing edges in the WT network, the knockout leaves the network unchanged and a warning says the results are only numerical noise. Target genes without outgoing edges are also reported when only some of them lack edges.
- **Transcriptome-wide perturbation** (`transcriptomeWide = TRUE`) — Builds the WT network **once**, then knocks out every gene in the network (or a user-supplied subset passed through `gKO`). It returns a matrix of perturbation distances, and a matrix of predicted directions, instead of a single differential regulation table. By default it uses the heat manifold alignment, which computes the heat kernel of the WT network once and reads every knockout from it; set `ma_method = "manifold"` to run one manifold alignment per perturbed gene instead, whose running time scales with the number of genes.

The comparison between the WT and KO networks is selected with `ma_method`. Both methods are available in both modes: `"manifold"` (the default for single and multi-gene knockouts) runs the non-linear manifold alignment of the WT and KO networks, and `"heat"` (the default for transcriptome-wide perturbation) uses the heat manifold alignment. The heat manifold alignment ranks the perturbed genes similarly to the manifold alignment (in the benchmark, AUROC for detecting the genes that change of 0.557 vs 0.555 in one atlas and 0.562 vs 0.580 in another), without recomputing an alignment for each knockout; use it when many knockouts are needed, and the manifold alignment when the best ranking for a few knockouts matters.

## Direction of the Response

When `dr_direction = TRUE` (the default), `scTenifoldKnk()` predicts whether each gene goes **up** or **down** after the knockout, using only the WT expression data (`knockoutDirection()`): the heat kernel of the gene-gene correlation matrix of log1p(CPM) expression is diffused from the knocked-out gene(s), and the sign of the result is the predicted direction.

Benchmarked against bulk knockdown/knockout profiles of the same cell lines, the predicted direction is correct more often than chance in two of the three single-cell atlases tested, and it is stable across random seeds and the number of genes used (direction AUROC 0.56 on 892 knockout experiments in the largest atlas; there, the genes with the largest `directionScore` move in the predicted direction 71% of the time, compared with 55% expected by chance). Three limitations apply:

- The predicted direction mostly reflects the response shared by most knockdowns along the dominant WT expression program, rather than regulation specific to the knocked-out gene. It is weakest for transcription factor knockouts.
- Its accuracy varies between cell types (strong in some cell lines, close to chance in others).
- It depends on the WT data set: with the third atlas, the same knockouts were predicted at close to chance level (direction AUROC 0.51), and the genes that change were also detected less accurately.

The magnitude (which genes respond) and the direction are reported separately (`distance`/`p.value` and `direction`/`directionScore`), so the direction can be used or ignored independently of the differential regulation statistics.

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
| `ma_method` | `NULL` | `"manifold"` (non-linear manifold alignment for each knockout) or `"heat"` (heat manifold alignment, one heat kernel for all knockouts). `NULL` uses `"manifold"` for single and multi-gene knockouts and `"heat"` for transcriptome-wide perturbation. |
| `ma_heatT` | `10` | Diffusion time of the heat kernel of the WT network (`ma_method = "heat"`). |
| `dr_empiricalNull` | `FALSE` | Assign differential regulation p-values using Efron's empirical null (via `locfdr`) instead of the theoretical chi-square null. |
| `dr_direction` | `TRUE` | Predict the direction (up/down) of the change of each gene with `knockoutDirection()`. |
| `dr_directionT` | `5` | Diffusion time of the correlation heat kernel used to predict the direction. |
| `nCores` | `parallel::detectCores()` | Number of cores used by the manifold alignment. |
| `seed` | `1` | Seed set before each random stage, so results are reproducible and do not depend on the caller's RNG state (which is restored on exit). Use different values to assess run-to-run variability, or `NULL` to use the caller's RNG state (e.g. a previous `set.seed()`). |

### `scQC()`

Standalone single-cell quality control. Arguments: `X` (raw counts matrix), `minLibSize`, `removeOutlierCells`, `minPCT`, `maxMTratio`, and an optional `label` for progress messages. Returns a `dgCMatrix` containing the cells and genes that pass the filters.

### `dRegulation()`

Differential regulation testing from a manifold alignment. Arguments: `manifoldOutput` (the labeled `manifoldAlignment` matrix, `X_` genes followed by `Y_` genes in the same order), `gKO` (the knocked-out genes, left out of the expectation used for the fold-changes; `scTenifoldKnk()` passes them) and `empiricalNull` (if `TRUE`, estimate the null distribution of the Z-scores with Efron's empirical null via the `locfdr` package instead of the theoretical chi-square null). An optional `direction` (named numeric vector of direction scores, e.g. a row of `knockoutDirection()`) adds the `direction` and `directionScore` columns. Returns the differential regulation table described under [Output](#output).

### `heatKernel()`

Spectral heat kernel of a gene-gene matrix, `H = sum_k exp(t (lambda_k / lambda_max - 1)) v_k v_k'`, from the eigendecomposition of its symmetric part. Arguments: `X` (square matrix), `t` (diffusion time; `t = 0` returns the identity) and `symmetric`.

### `hkManifoldAlignment()`

Heat manifold alignment of the WT network for one or many knockouts. The heat kernel of the WT network is computed once and the knockout of each gene `x` is read from it: `x` loses its mean log1p(CPM) expression and the change diffuses over the network, `delta_g = -mean(x) / sd(x) * H[x, g] * sd(g)`. Returns the matrix of perturbation distances `|delta_g|` (knockouts by genes). Arguments: `WT` (the WT network), `X` (WT raw counts), `gKO` (genes or a list of gene sets to knock out; `NULL` for every gene), `t` (default `10`) and an optional pre-computed kernel `H`.

### `knockoutDirection()`

Predicted direction of the response of each gene to a knockout, from the heat kernel of the WT gene-gene correlation matrix (see [Direction of the Response](#direction-of-the-response)). Arguments: `X` (WT raw counts), `gKO`, `genes` (genes to score) and `t` (default `5`). Returns a matrix of direction scores (positive: up, negative: down), knockouts by genes.

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

- **`tensorNetworks`** — Weight-averaged denoised gene regulatory network after CP tensor decomposition, containing:
  - `WT`: The network for the wild-type condition (a `Matrix` object). The knocked-out network is the same network with the rows of `gKO` set to 0; it is not returned, which roughly halves the size of the output (rebuild it with `KO <- output$tensorNetworks$WT; KO[gKO, ] <- 0`).
- **`manifoldAlignment`** — A data frame of low-dimensional features from the non-linear manifold alignment, with 2 × *n* genes rows and *d* columns (default *d* = 2). Only returned when `ma_method = "manifold"`.
- **`diffRegulation`** — A data frame with eight columns (six when `dr_direction = FALSE`):
  - `gene`: Gene identifier.
  - `distance`: Euclidean distance between the gene's coordinates in the two conditions.
  - `Z`: Z-score after Box-Cox power transformation.
  - `FC`: Fold change of the squared distance with respect to the expectation, the mean squared distance of the genes that were not knocked out.
  - `p.value`: P-value from the chi-square distribution with one degree of freedom, or from Efron's empirical null when `dr_empiricalNull = TRUE`.
  - `p.adj`: Adjusted p-value (Benjamini & Hochberg FDR correction).
  - `direction`: Predicted direction of the change of the gene, `"up"` or `"down"` (the knocked-out genes are `"down"` by construction).
  - `directionScore`: Signed direction score from `knockoutDirection()` (`NA` for the knocked-out genes).

### Transcriptome-wide mode (`transcriptomeWide = TRUE`)

`scTenifoldKnk()` returns a list with three elements:

- **`tensorNetworks`** — A list with the WT weight-averaged denoised gene regulatory network (`WT`).
- **`perturbationDistances`** — A numeric matrix of perturbation distances (heat manifold alignment by default, manifold-alignment distances with `ma_method = "manifold"`). Rows are the perturbed genes, columns are all genes in the WT network, and each entry is the distance of a gene under the corresponding knockout.
- **`perturbationDirections`** — A numeric matrix with the same dimensions and the predicted direction of each gene under each knockout (`1` up, `-1` down, `0` undetermined). Only returned when `dr_direction = TRUE`.

## Running Time

Running time grows mainly with the number of genes. The number of cells matters little, because each network is built from a fixed-size subsample of cells (`nc_nCells`). Single knockout benchmarks with the default parameters (10 networks of 500 cells) on simulated counts, measured on an Apple M2 Pro (16 GB RAM) with R 4.5 and its reference BLAS:

| Cells | Genes | Time |
|------:|------:|-----:|
| 300 | 1,000 | 17 s |
| 1,000 | 1,000 | 17 s |
| 1,000 | 5,000 | 4.0 min |
| 2,500 | 5,000 | 3.5 min |
| 5,000 | 5,000 | 3.7 min |

Before version 1.1.1 (scTenifoldNet 1.4.1), networks were built by fitting one SVD per gene, and the earlier benchmarks in this README (measured on a different machine) reported about 3 hours for 5,000 genes.

### Memory

Peak memory grows with the square of the number of genes, because the networks are stacked into a tensor of genes x genes x `nc_nNet` entries. Estimated peak memory with 10 networks:

| Genes | Peak memory |
|------:|------------:|
| 1,000 | 0.5 GB |
| 2,000 | 1.5 GB |
| 5,000 | 8.6 GB |
| 10,000 | 34 GB |
| 15,000 | 76 GB |

`scTenifoldKnk()` compares this estimate with the memory available after quality control and warns, before building the networks, when it does not fit, reporting the largest number of genes that does. The estimate can also be checked beforehand with `scTenifoldNet::checkMemory(nGenes, nNet = 10, nConditions = 1)`.

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

# Accessing the WT gene regulatory network, and rebuilding the KO network
dim(output$tensorNetworks$WT)
KO <- output$tensorNetworks$WT
KO['ng10', ] <- 0

# Accessing the manifold alignment result
head(output$manifoldAlignment)

# Differential regulation results — top perturbed genes, with the predicted direction
head(output$diffRegulation, n = 10)

# The same knockout with the heat manifold alignment
heatOutput <- scTenifoldKnk(X, gKO = "ng10", ma_method = "heat", qc_minLibSize = 30)
head(heatOutput$diffRegulation, n = 10)

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

Knock out every gene in the WT network and collect the perturbation distances and directions (heat manifold alignment by default):

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

# Predicted directions (1 up, -1 down)
twOutput$perturbationDirections[1:5, 1:5]

# The manifold alignment for every perturbed gene remains available
twManifold <- scTenifoldKnk(X, transcriptomeWide = TRUE, ma_method = "manifold", qc_minLibSize = 30)
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
