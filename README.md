
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toxassay <a href="https://dplyr.tidyverse.org"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dplyr)](https://cran.r-project.org/package=dplyr)
[![R-CMD-check](https://github.com/tidyverse/dplyr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tidyverse/dplyr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/tidyverse/dplyr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tidyverse/dplyr?branch=main)
<!-- badges: end -->

## Overview

ToxAssay is an R-based software package designed for the comprehensive
evaluation of drug-induced toxicity utilizing complex toxicogenomics
databases. The package offers a comprehensive suite of functions,
encompassing:

- Automated acquisition, preprocessing, and annotation of data by
  compound name for efficient management.
- Identification of molecular markers for targeted toxicity, such as
  differentially expressed genes (DEGs),adverse outcome pathways (AOPs),
  functional pathways, and Protein-Protein Interaction (PPI) networks,
  etc.
- Development of optimized machine-learning classifiers for predicting
  the targeted toxicity in test samples of compounds.

Workflow of ToxAssay in five main steps:
<img src="man/figures/workflow.png" width="890"/>

## Installation

``` r
# Install released version from CRAN
install.packages("toxassay")
```

To get a bug fix or to use a feature from the development version, you
can install the development version of toxassay from GitHub.

``` r
# install.packages("pak")
pak::pak("ringku09/toxassay")
```

## Usages

### Download and preprocessing analysis ready data

In the current version of **ToxAssay**, users can download and process
raw compound perturbation data from the [Open
TG-GATEs](https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html)
and [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) databases,
as well as curated relational data from the [CTD](https://ctdbase.org).

Available compounds in databases.

``` r
# Compounds in Open TG-GATEs database
# Data available for compounds in Rat in-vivo liver Single data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")

# Compounds in Open DrugMatrix database
# Data available for compounds in `Liver`
drugmatrix_compounds(tissue = "Liver")
```

``` r
# Get data from Open TG-GATEs
get_tggates(compounds = c("2-nitrofluorene", "N-nitrosomorpholine"), species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")

# Get data from Open DrugMatrix
get_drugmatrix(compounds = c("Ethanol", "Estriol"), tissue = "Liver")

# Get data from CTD
get_ctd(compounds = c("Ethanol", "Estriol"), gene_parm = "cgixns", disease_parm = "diseases")
```

<img src="man/figures/downloading.png" width="890"/>

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/ringku09/toxassay/issues). For questions and
other discussion, please send an [E-mail](mailto:ringku_740@yahoo.com).
