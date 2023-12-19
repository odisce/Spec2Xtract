
# Spec2Xtract <img src="man/figures/logo.svg" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/odisce/Spec2Xtract/branch/main/graph/badge.svg)](https://app.codecov.io/gh/odisce/Spec2Xtract?branch=main)
[![R-CMD-check](https://github.com/odisce/Spec2Xtract/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/odisce/Spec2Xtract/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Spec2Xtract** is an [R](https://www.r-project.org/) package to extract MS2 spectra from 
[.mzML](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#mzML) files or directly from .raw files 
acquired on [Thermo](https://www.thermofisher.com) instrument by using [mzR](https://github.com/sneumann/mzR/) 
or [rawrr](https://github.com/fgcz/rawrr) respectively. The package is compatible with multi-energy acquisition, 
extracting one spectra for every energy used in the file. Spectra annotation can be performed by using 
[**Spec2Annot**](https://github.com/odisce/Spec2Annot). An implementation on 
[Galaxy](https://workflow4metabolomics.usegalaxy.fr/) of those functionalities is under progress.

## Installation

You can install the development version of Spec2Xtract from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("odisce/Spec2Xtract")
```

## Example

### Extract MS2 spectra from a .raw file

In progress

