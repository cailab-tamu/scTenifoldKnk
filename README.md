scTenifoldKnk
=============

Install:
-------
This package is under active development, you can install **scTenifoldKnk**, using the following command:

```{R}
library(remotes)
install_github('cailab-tamu/scTenifoldKnk')
library(scTenifoldKnk)
```

Available functions:
--------------------

|Code| Function |
|:-|:-|
|scTenifoldKnk||

Input:
--------
The required input for **scTenifoldKnk** is an expression matrix with genes in the rows and cells (barcodes) in the columns. Data is expected to be _not normalized_.

Running time:
--------
The running time of scTenifoldKnk is largely dependent on how long it takes to construct scGRNs from subsampled expression matrices. Time increases proportional to the number of cells and genes in the dataset used as input. Below is a table of running times under different scenarios:

| Number of Cells | Number of Genes | Running Time |
|-----------------|-----------------|--------------|
| 300             | 1000            | 3.45 min     |
| 1000            | 1000            | 4.25 min     |
| 1000            | 5000            | 171.88 min (2 h 51.6 min) |
| 2500            | 5000            | 175.29 min (2 h 55.3 min) |
| 5000            | 5000            | 188.88 min (3 h 8.9 min) |
| 5000            | 7500            | 189.51 min (3 h 9.5 min)  |
| 7500            | 5000            | 615.45 min (10 h 15.5 min) |
| 7500            | 7500            | 616.12 min (10 h 16.1 min)  |


Output:
--------
The output of **scTenifoldKnk** is a list with 3 slots as follows: 
  * **tensorNetworks**: The computed weight-averaged denoised gene regulatory networks after CANDECOMP/PARAFAC (CP) tensor decomposition. It includes two slots with:
    * **X**: The constructed network for the _X_ sample.
    * **Y**: The constructed network for the _Y_ sample.
  * **manifoldAlignment**: The generated low-dimensional features result of the non-linear manifold alignment. It is a data frame with _2 times the number of genes_ in the rows and _d_ (default= 2) dimensions in the columns
  * **diffRegulation**: The results of the differential regulation analysis. It is a data frame with 6 columns as follows:
    * **gene**: A character vector with the gene id identified from the manifoldAlignment output.
    * **distance**: A numeric vector of the Euclidean distance computed between the coordinates of the same gene in both conditions.
    * **Z**: A numeric vector of the Z-scores computed after Box-Cox power transformation.
    * **FC**: A numeric vector of the FC computed with respect to the expectation.
    * **p.value**: A numeric vector of the p-values associated to the fold-changes, probabilities are asigned as P[X > x] using the Chi-square distribution with one degree of freedom.
    * **p.adj**: A numeric vector of adjusted p-values using Benjamini & Hochberg (1995) FDR correction.
