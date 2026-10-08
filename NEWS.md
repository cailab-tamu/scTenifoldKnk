# scTenifoldKnk 2.0.0

## New features

- **Direction of the response.** `knockoutDirection()` predicts whether each gene goes up or down after a knockout, from the heat kernel of the WT gene-gene correlation matrix of log1p(CPM) expression. `scTenifoldKnk()` reports it by default (`dr_direction = TRUE`, `dr_directionT = 5`) as the `direction` and `directionScore` columns of `diffRegulation`, and as `perturbationDirections` in transcriptome-wide mode. `dr_directionRegressLibSize = TRUE` (off by default) regresses log library size out of each gene first, for data sets in which nearly all genes are predicted down.
- **Heat manifold alignment.** `hkManifoldAlignment()` computes the heat kernel of the WT network once (`heatKernel()`) and reads every knockout from it, so the cost of a transcriptome-wide screen no longer scales with the number of knockouts. Selected with `ma_method = "heat"` (`ma_heatT = 10`).
- **`perturbationMap()`** ranks the knockouts of a transcriptome-wide run by the similarity of their signed profiles to a reference signature, such as a disease-vs-control comparison of the same cell type, and draws them on a map.

## Changes in behaviour

- `diffRegulation` has eight columns instead of six (`direction`, `directionScore`). Set `dr_direction = FALSE` for the 1.x table.
- Transcriptome-wide mode (`transcriptomeWide = TRUE`) uses the heat manifold alignment by default, so its distances differ from 1.x. Set `ma_method = "manifold"` to recover them. Single and multi-gene knockouts still use the manifold alignment by default and give the same distances, Z-scores and p-values as 1.1.5.

## Documentation

- README rewritten: quality-control recommendations (removal of ribosomal and mitochondrial genes), use in tissues and disease data, and the limitations of the predicted direction.

# scTenifoldKnk 1.1.5

- The knocked-out genes are left out of the expectation of the differential regulation test (#45).

# scTenifoldKnk 1.1.4

- Distances at the level of floating-point noise are ignored by `dRegulation()`.

# scTenifoldKnk 1.1.3

- Warning when the networks of the requested genes do not fit in memory (#10).

# scTenifoldKnk 1.1.2

- Multi-gene knockout: several genes in `gKO` are knocked out together (#47).
- `plotKO()` enrichment works without attaching enrichR (#43); data.frame input is rejected (#44).

# scTenifoldKnk 1.1.1

- `seed` argument for reproducible results (#48); warning when the knocked-out genes have no outgoing edges.

# scTenifoldKnk 1.1.0

- Transcriptome-wide perturbation mode and Efron's empirical null for the differential regulation test.
