# scTenifoldKnk plotKO() — Frequently Asked Questions

Common questions about the `plotKO()` visualization, edge-weight interpretation,
and the scope of differential-regulation inference.
Answers provided by the original developers.

---

## Q1: What do the red and blue edge colors in the `plotKO()` egocentric network plot represent?

Edge colors in the `plotKO()` visualization encode the **sign of the co-expression relationship** between the knocked-out gene and each of its neighbors in the wild-type (WT) single-cell gene regulatory network (scGRN):

```r
edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue')
```

Here, `W` is the edge weight drawn from the WT adjacency matrix (`X$tensorNetworks$WT`). The weight follows the **sign of the correlation** between the two genes:

- **Red edge (W > 0)** — positive co-expression between the KO gene and its neighbor in the WT scGRN.
- **Blue edge (W < 0)** — negative co-expression between the KO gene and its neighbor in the WT scGRN.

This interpretation is consistent with the PC regression framework used to construct the tensor network, where edge weights reflect the sign of the correlation between genes (see Figure 4 in [Cai et al. (2020), *Cells*](https://www.mdpi.com/2073-4409/9/1/14)).

---

## Q2: Does `plotKO()` indicate whether a perturbed gene will be *upregulated* or *downregulated* following the virtual knockout?

**No.** scTenifoldKnk identifies *which* genes are significantly differentially regulated (DR) after the virtual knockout, but it does not predict the **direction** of expression change (up vs. down) for those genes.

The fold-change (`FC`) values reported in the `diffRegulation` output are derived from a chi-square test and are therefore always positive — they represent the *magnitude* of perturbation, not its direction.

The reason directional prediction is not currently possible is that **a gene typically has more than one regulator**. Because scTenifoldKnk perturbs only one gene at a time, there is no way to infer the net expected change in expression of a downstream gene whose response depends on the balance of multiple regulators.

> **Empirical note:** The original scTenifoldKnk paper (Patterns, 2022) reports that significantly DR genes tend to be *downregulated* in real knockout experiments. This is an empirical observation, not a direct output of the tool, and may reflect the general effect of removing a positive regulator.

---

## Q3: Does the edge color in `plotKO()` provide any probabilistic hint about the direction of downstream regulation after knockout?

The correlation sign captured by edge color does carry **directional information about the regulatory relationship** between two genes in the WT network. In principle, a red edge (positive co-expression) between a KO gene and a downstream neighbor is *consistent with* the neighbor tending to be downregulated when the KO gene is removed — because the positive regulator has been eliminated.

However, this is **not a reliable prediction** at the level of individual genes. Because most genes are regulated by multiple inputs, the actual direction of expression change after a knockout depends on the combined effect of all regulators, not just the one being perturbed.

Providing a proper directional perturbation estimate — for example, by modeling network-based diffusion where the in-degree of each gene is used to standardize information flow and a signal is propagated until the network reaches steady state — is a planned future enhancement to scTenifoldKnk.

> **Summary:** Edge color gives directional context about the *WT regulatory relationship*, but should not be used alone to predict up- or downregulation after knockout. Use it alongside the `diffRegulation` output and published experimental validation data.

---

## References

- Osorio et al. (2022). *scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation.* [Patterns.](https://doi.org/10.1016/j.patter.2022.100434)
- Cai et al. (2020). *Single-Cell Expression Variability Implies Cell Function.* [Cells.](https://www.mdpi.com/2073-4409/9/1/14)
