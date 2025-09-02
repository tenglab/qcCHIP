# qcCHIP

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17037970.svg)](https://doi.org/10.5281/zenodo.17037970)

A package to perform quality filtering of clonal hematopoiesis (CH) mutations using cohort-specific characteristics. Four types of filtering metrics are implemented: technical metrics, functional metrics, individual metrics and population metrics. It also allows users to determine the optimal values of the following numerical metrics based on permutation analysis: VAF, DP, mutation prevalence, SOR, and SAF/SAR. 

## Installation

qcCHIP is an R package with its source code documented in this [repository](https://github.com/tenglab/qcCHIP.git).


```R
# Install qcCHIP
library(devtools)
devtools::install_github("https://github.com/tenglab/qcCHIP.git")
```

Or, Install with vignettes and dependencies.

```R
# Install qcCHIP with vignettes and dependencies
devtools::install_github("https://github.com/tenglab/qcCHIP.git",build_vignettes = TRUE)
```

## Using qcCHIP
First load qcCHIP,
```R
library(qcCHIP)
```
Then, follow the [**User's Guide**](https://github.com/tenglab/qcCHIP/blob/main/vignettes/qcCHIP.pdf)
to load qcCHIP. In the guide, we detail how to use qcCHIP functions to determine values of filtering metrics and to filter CH
mutations with the four types of metrics.

The users are also encouraged to refer to the help pages of R functions in this package. 

## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/qcCHIP/issues).
