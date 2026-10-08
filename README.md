# scTenifoldKnk

[![R-CMD-check](https://github.com/cailab-tamu/scTenifoldKnk/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cailab-tamu/scTenifoldKnk/actions/workflows/R-CMD-check.yaml)
[![CRAN](https://www.r-pkg.org/badges/version/scTenifoldKnk)](https://CRAN.R-project.org/package=scTenifoldKnk)
[![License: GPL (>=2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D2%29-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

**scTenifoldKnk** performs virtual knockout experiments from single-cell RNA-seq data of wild-type (WT) cells. It builds a gene regulatory network (GRN) from the WT cells, removes the outgoing edges of the target gene, and compares the knocked-out network with the WT network to identify the genes whose regulation changes. Version 2.0 also predicts the **direction** (up or down) of each gene's response and scores **every gene of the network as a knockout** from a single heat kernel.

Other implementations: [scTenifoldpy](https://github.com/qwerty239qwe/scTenifoldpy) (Python, same results) and [scGEAToolbox](https://github.com/jamesjcai/scGEAToolbox) (MATLAB).

## Installation

```r
install.packages("scTenifoldKnk")                       # CRAN
remotes::install_github("cailab-tamu/scTenifoldKnk")    # development version
```

## What is new in 2.0

| Feature | Function | Argument of `scTenifoldKnk()` |
|:--|:--|:--|
| Direction of the response of each gene | `knockoutDirection()` | `dr_direction = TRUE` (default) |
| Heat manifold alignment: all knockouts from one heat kernel of the WT network | `heatKernel()`, `hkManifoldAlignment()` | `ma_method = "heat"` |
| Comparison of transcriptome-wide knockouts with a reference signature | `perturbationMap()` | |

Changes in behaviour: `diffRegulation` gains the `direction` and `directionScore` columns, and transcriptome-wide mode uses the heat manifold alignment by default. `dr_direction = FALSE` and `ma_method = "manifold"` restore the 1.x output. See [NEWS.md](NEWS.md).

## Method

| Step | Function | Description |
|:--:|:--|:--|
| 1 | `scQC` | Cell and gene filters: library size, outliers, detection rate, mitochondrial fraction |
| 2 | `cpmNormalization` | Counts-per-million normalization |
| 3 | `makeNetworks` | Principal component regression networks (`pcNet`) on subsamples of cells |
| 4 | `tensorDecomposition` | CANDECOMP/PARAFAC decomposition of the stacked networks (denoising) |
| 5 | `strictDirection` | Keeps the stronger direction of each pair of edges |
| 6 | `manifoldAlignment` / `hkManifoldAlignment` | Comparison of the WT and KO networks |
| 7 | `dRegulation` | Box-Cox transformed distances tested against a chi-square null, with the predicted direction |

**Knockout.** The rows of the knocked-out genes in the WT network are set to zero.

**Manifold alignment** (`ma_method = "manifold"`, default for single and multi-gene knockouts) embeds the WT and KO networks in a shared low-dimensional space; the perturbation of each gene is the distance between its two embeddings. One alignment is computed per knockout.

**Heat manifold alignment** (`ma_method = "heat"`, default for transcriptome-wide mode) computes the heat kernel of the WT network once, `H = sum_k exp(t (lambda_k / lambda_max - 1)) v_k v_k'` (`ma_heatT = 10`). Knocking out gene *x* removes its mean log1p(CPM) expression, which diffuses over the network: `delta_g = -mean(x) / sd(x) * H[x, g] * sd(g)`. The perturbation distance of gene *g* is `|delta_g|`. It ranks genes similarly to the manifold alignment, and its cost does not grow with the number of knockouts.

**Direction** (`knockoutDirection()`) uses the same diffusion on the heat kernel of the WT gene-gene correlation matrix of log1p(CPM) expression (`dr_directionT = 5`); the sign of the result is the predicted direction. The magnitude (`distance`, `p.value`) and the direction (`direction`, `directionScore`) are reported separately and can be used independently.

## Usage

The input is a raw count matrix (`matrix` or `dgCMatrix`), genes in rows and cells in columns.

```r
library(scTenifoldKnk)

set.seed(1)
X <- matrix(rnbinom(100 * 2000, size = 20, prob = 0.98), ncol = 2000)
rownames(X) <- c(paste0("ng", 1:90), paste0("mt-", 1:10))

# Single knockout
out <- scTenifoldKnk(X, gKO = "ng10", qc_minLibSize = 30)
head(out$diffRegulation)          # gene, distance, Z, FC, p.value, p.adj, direction, directionScore
plotKO(out, gKO = "ng10")

# Multi-gene knockout (genes knocked out together)
dko <- scTenifoldKnk(X, gKO = c("ng10", "ng20"), qc_minLibSize = 30)

# Transcriptome-wide: every gene knocked out separately
tw <- scTenifoldKnk(X, transcriptomeWide = TRUE, qc_minLibSize = 30)
tw$perturbationDistances[1:5, 1:5]   # knockouts x genes
tw$perturbationDirections[1:5, 1:5]  # 1 up, -1 down

# Knockouts ranked by similarity to a reference signature (e.g. disease vs control log fold changes)
signature <- setNames(rnorm(ncol(tw$perturbationDistances)), colnames(tw$perturbationDistances))
pm <- perturbationMap(tw, signature = signature, genes = c("ng10", "ng20"))
```

Every argument is documented in `?scTenifoldKnk`. The network size (`nc_nNet = 10` networks of `nc_nCells = 500` cells, `nc_q = 0.9`) and the tensor rank (`td_K = 3`) follow the original publication.

## Output

**Single or multi-gene knockout** (`transcriptomeWide = FALSE`):

- `tensorNetworks$WT`: denoised WT network. The KO network is the same matrix with the rows of `gKO` set to zero.
- `manifoldAlignment`: coordinates of the WT (`X_`) and KO (`Y_`) genes (manifold alignment only).
- `diffRegulation`: one row per gene, sorted by significance.

  | Column | Description |
  |:--|:--|
  | `distance` | Perturbation distance between the WT and KO networks |
  | `Z` | Z-score of the Box-Cox transformed distance |
  | `FC` | Squared distance relative to the mean squared distance of the genes not knocked out |
  | `p.value`, `p.adj` | Chi-square (df = 1) p-value, or Efron's empirical null with `dr_empiricalNull = TRUE`; Benjamini-Hochberg adjustment |
  | `direction`, `directionScore` | Predicted change (`"up"`/`"down"`) and its signed score |

**Transcriptome-wide** (`transcriptomeWide = TRUE`): `tensorNetworks$WT`, `perturbationDistances` and `perturbationDirections` (knockouts x genes). Pass a subset of genes in `gKO` to perturb only those, each one separately.

## Recommendations

**Quality control.** Remove ribosomal and mitochondrial genes before building the networks. Their dense co-expression modules otherwise dominate the network and the top of every knockout. Apply the mitochondrial read filter first (`scQC()`), then:

```r
riboMito <- grepl("^(RP[LS][0-9]|RPLP[0-9]|RPSA|MRP[LS][0-9]|MT-|MTND[0-9]|MTCO[0-9]|MTATP[0-9]|MTCYB|MTRNR2L)",
                  toupper(rownames(X)))
X <- X[!riboMito, ]
```

**Tissues and disease data.**

- Build one network per cell type or state (about 500 cells or more), so that cell identity does not dominate the co-expression.
- Correct ambient RNA (e.g. DecontX) when other cell types contribute contaminating transcripts.
- To prioritize candidate genes, run the transcriptome-wide mode on the cell type of interest and compare each knockout with a disease-vs-control signature of the same cell type (`perturbationMap()`). A known loss-of-function causal gene, when available, is a positive control: its knockout in healthy cells should reproduce the patient signature.
- Genes perturbed by most knockouts reflect a shared response rather than the knocked-out gene; the transcriptome-wide distance matrix provides this background.

**Seeds.** Network construction subsamples cells. Results are reproducible for a given `seed` (default `1`) and independent of the caller's RNG state; averaging two or three seeds smooths the ranking of a single knockout. `seed = NULL` uses the caller's RNG.

## Limitations

- The predicted direction mainly reflects the response shared by most perturbations along the dominant WT expression program, rather than regulation specific to the knocked-out gene. Its accuracy varies between cell types and data sets.
- In single-cell data, log1p(CPM) expression can still follow sequencing depth, which makes nearly all genes correlate positively and predicts nearly every gene down. `dr_directionRegressLibSize = TRUE` removes that axis; it is off by default because the depth-associated component can be biological in cell lines.
- `perturbationMap()` scores similarity by cosine, which is sensitive to the overall sign of a profile. Check the gene-level agreement of top-ranked knockouts before interpreting them.
- Virtual knockouts prioritize candidates for experimental testing; they do not replace it.

## Running time and memory

Running time grows with the number of genes, not cells, because each network uses a fixed-size subsample of cells. Peak memory grows with the square of the number of genes (genes x genes x `nc_nNet` tensor); `scTenifoldKnk()` warns when the estimate exceeds the available memory (`scTenifoldNet::checkMemory()`).

| Genes | Single knockout | Transcriptome-wide, heat | Transcriptome-wide, manifold | Peak memory |
|--:|--:|--:|--:|--:|
| 1,000 | 17 s | 13 s | ~3 min | 0.5 GB |
| 3,000 | | 2.4 min | ~39 min | |
| 5,000 | 3.7 min | 8.2 min | ~3.6 h | 8.6 GB |
| 10,000 | | | | 34 GB |

Single knockouts: 1,000 to 5,000 cells, Apple M2 Pro, R 4.5 reference BLAS. Transcriptome-wide: network, heat kernel and directions for all genes on an Apple M4; manifold times are extrapolated from the measured time per knockout.

## Citation

Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Tian Y, Chapkin RS, Huang JZ, Cai JJ. scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation. *Patterns* 3(3):100434 (2022). [doi:10.1016/j.patter.2022.100434](https://doi.org/10.1016/j.patter.2022.100434)

```bibtex
@article{osorio2022sctenifoldknk,
  title   = {scTenifoldKnk: An Efficient Virtual Knockout Tool for Gene Function
             Predictions via Single-Cell Gene Regulatory Network Perturbation},
  author  = {Osorio, Daniel and Zhong, Yan and Li, Guanxun and Xu, Qian and Yang, Yongjian and
             Tian, Yanan and Chapkin, Robert and Huang, Jianhua Z. and Cai, James J.},
  journal = {Patterns}, year = {2022}, volume = {3}, number = {3}, pages = {100434},
  doi     = {10.1016/j.patter.2022.100434}
}
```

---

&copy; The Texas A&M University System.
