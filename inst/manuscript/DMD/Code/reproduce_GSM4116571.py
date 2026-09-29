"""Reproduce the Dmd virtual knockout of the manuscript (GSM4116571) from the raw data.

Results/scTenifoldKnk_SRS4245406.csv was computed in April 2020 with the
scTenifoldNet code of that time (commit acc68ed) on R < 3.6. The current
packages compute the same PC networks, but four details of that run cannot
be set through their parameters, so this script sets them:

1. Cells: R < 3.6 sampled with the "Rounding" sampler, floor(n * runif()),
   so set.seed(1) picked other cells than it does now.
2. Tensor decomposition: the tensor had 4 modes (genes x genes x 1 x
   networks). R drew one more random 1 x K factor before the network factor.
   After its first update that factor is a sign per component, so the CP-ALS
   is the 3-mode one started from the network factor times that draw.
3. Raw counts (no CPM), K = 5, rounded to 1 decimal, the Dmd row set to 0
   (no transposition), and d = 2 in the manifold alignment.
4. Ribosomal and mitochondrial genes were removed after the manifold
   alignment, not before as in the commented code of DMD_DataProcessing.R.

Last, the knocked-out gene is left out of the expectation of the chi-square
test, as scTenifoldKnk did until 1.0.3 and again from 1.1.5 (scTenifoldpy
0.5.1).

The WT network is identical to the one in Results/GSM4116571.RData, and the
190 genes with p.adj < 0.05 are the same as in the published table, with the
same top genes. Only the distances below ~1e-11, far below the floating-point
noise level, differ; through the Box-Cox power they change the Z-scores by up
to about 0.3.

Requires scTenifoldpy >= 0.5.1 and about 8 GB of memory; it takes 3 to 5
minutes on a 12-core laptop.
Run it from this directory:

    python reproduce_GSM4116571.py
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

import scTenifold.core._decomposition as decomposition
import scTenifold.core._networks as networks
from scTenifold.core._decomposition import tensor_decomp
from scTenifold.core._networks import d_regulation, make_networks, manifold_alignment
from scTenifold.core._rng import RRandom

DATA = Path("../Data")
RESULTS = Path("../Results")
K = 5


def read_data():
    X = scipy.io.mmread(DATA / "GSM4116571_Qc_matrix.mtx.gz").tocsc()
    genes = pd.read_csv(DATA / "GSM4116571_Qc_features.tsv.gz", sep="\t", header=None)[1].to_numpy()
    with gzip.open(DATA / "GSM4116571_Qc_barcodes.tsv.gz", "rt") as f:
        cells = f.read().split()
    return X, genes, np.array(cells)


def prediction_interval(x, y, level=0.95):
    """``predict(lm(y ~ x), interval = 'prediction')`` in R."""
    n = len(x)
    slope, intercept = np.polyfit(x, y, 1)
    fit = intercept + slope * x
    sigma = np.sqrt(np.sum((y - fit) ** 2) / (n - 2))
    se = sigma * np.sqrt(1 + 1 / n + (x - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    t = stats.t.ppf((1 + level) / 2, n - 2)
    return fit - t * se, fit + t * se


def sc_qc_2020(X, genes, cells, mt_threshold=0.1, min_lib_size=1000):
    """scQC of github.com/dosorio/utilities (singleCell/scQC.R) used in 2020."""
    lib_size = np.asarray(X.sum(axis=0)).ravel()
    X, cells = X[:, lib_size >= min_lib_size], cells[lib_size >= min_lib_size]
    lib_size = np.asarray(X.sum(axis=0)).ravel().astype(float)
    n_genes = np.asarray((X != 0).sum(axis=0)).ravel().astype(float)
    is_mt = pd.Series(genes).str.upper().str.match("^MT-").to_numpy()
    mt_counts = np.asarray(X[is_mt].sum(axis=0)).ravel().astype(float)
    g_lwr, g_upr = prediction_interval(lib_size, n_genes)
    m_lwr, m_upr = prediction_interval(lib_size, mt_counts)
    keep = ((mt_counts > m_lwr) & (mt_counts < m_upr) & (n_genes > g_lwr) & (n_genes < g_upr)
            & (mt_counts / lib_size <= mt_threshold) & (lib_size < 2 * lib_size.mean()))
    return X[:, keep], cells[keep]


class RRandom2020(RRandom):
    """R's random numbers as in the April 2020 run."""

    def __init__(self, seed=None):
        super().__init__(seed)
        self._rnorm_calls = 0

    def _unif_index(self, n):
        # sample.kind = "Rounding", the default before R 3.6.0
        return int(np.floor(n * self.unif_rand(1)[0]))

    def rnorm(self, n):
        # cp_decomposition draws the gene, gene and network factors; the 4-mode
        # tensor of 2020 drew a 1 x K factor before the network factor
        self._rnorm_calls += 1
        if self._rnorm_calls == 3:
            c = super().rnorm(K)
            D = super().rnorm(n).reshape((n // K, K), order="F") * c
            return D.reshape(-1, order="F")
        return super().rnorm(n)


def main():
    networks.RRandom = RRandom2020
    decomposition.RRandom = RRandom2020

    X, genes, cells = read_data()
    X, cells = sc_qc_2020(X, genes, cells)
    # rowMeans(X != 0) > 0.05; divided as in R, so genes in exactly 5% of the cells are left out
    expressed = np.asarray((X != 0).sum(axis=1)).ravel() / X.shape[1] > 0.05
    X, genes = X[expressed], genes[expressed]
    print(f"After QC: {X.shape[1]} cells, {X.shape[0]} genes")
    counts = pd.DataFrame(X.toarray(), index=genes, columns=cells)

    # scTenifoldNet::makeNetworks(countMatrix, q = 0.9) after set.seed(1)
    nets = make_networks(counts, n_nets=10, n_samp_cells=500, n_comp=3, q=0.9, random_state=1,
                         backend="joblib-threading", n_jobs=8)
    # scTenifoldNet::tensorDecomposition(WT) with the defaults of 2020
    WT = tensor_decomp([n.toarray() for n in nets], list(genes), K=K, n_decimal=1, tol=1e-5,
                       max_iter=1000, random_state=1)
    del nets
    KO = WT.copy()
    KO.loc["Dmd", :] = 0
    MA = manifold_alignment(WT, KO, d=2)

    keep = ~pd.Series(genes).str.contains(r"^Rp[0-9]+|^Rpl|^Rps|^Mt-", case=False).to_numpy()
    MA = MA[np.concatenate([keep, keep])]
    DR = d_regulation(MA, ko_genes=["Dmd"])
    DR.to_csv("scTenifoldKnk_GSM4116571_reproduced.csv", index=False)

    published = pd.read_csv(RESULTS / "scTenifoldKnk_SRS4245406.csv")
    sig = set(DR["Gene"][DR["adjusted p-value"] < 0.05])
    published_sig = set(published["gene"][published["p.adj"] < 0.05])
    print(f"Genes with p.adj < 0.05: {len(sig)} reproduced, {len(published_sig)} published, "
          f"{len(sig & published_sig)} in both")
    print("Top 10:", ", ".join(DR["Gene"][:10]))


if __name__ == "__main__":
    main()
