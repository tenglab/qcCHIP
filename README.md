# qcCHIP

[![DOI](https://zenodo.org/badge/960167508.svg)](https://doi.org/10.5281/zenodo.16421861)

A package of helping user select optimal VAF, DP, population metrics, SOR, and SAF/SAR based on a permutation approach. It also provide a function to extract CHIP candidates based on variety of metrics.

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
to load qcCHIP. In the guide, we detail how to query through qcCHIP
and to visualize SE signatures across cancers.

The users are also encouraged to refer to the help pages of R functions in this package. 

## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/qcCHIP/issues).
