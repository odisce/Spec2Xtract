
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

## Status

  [x] Main steps  
  [x] targets pipeline as a factory  
  [ ] Documentation  
  [ ] Tests  
  [ ] Vignette  

## Installation

You can install the development version of Spec2Xtract from [GitHub](https://github.com/) with:

| method | command |
|--|--|
| with `remotes` | `r remotes::install_github("odisce/Spec2Xtract")` |
| with `devtools` | `r devtools::install_github("odisce/Spec2Xtract")` |
| with `renv` | `r renv::install("github::odisce/Spec2Xtract")` |
| with `pak` | `r pak::pkg_install("odisce/Spec2Xtract")` |

`Spec2Xtract` use [`rawrr`](https://github.com/fgcz/rawrr) ro read *.raw* files directly. To be used it needs to 
be installed with the following commands:

```{r}
rawrr::installRawFileReaderDLLs()
rawrr::installRawrrExe()
```


## Example

`Spec2Xtract` provide a [`targets`](https://github.com/ropensci/targets) factory to run a full analysis pipeline
from a list of compound(s) and a list of path to *.raw* file(s) as shown below. [`targets`](https://github.com/ropensci/targets)
is a great package to manage analytical workflow and this implementation leverage all the functionalities it provides
(Caching, vizualisation, CPUs parallelization, HPC deployment, etc...).

First, create a new folder to store your project and open an R session inside: 

```{r}
## Load Spec2Xtract library
require("Spec2Xtract")

## Here we use a temporary directory
temp_wd <- tempdir()
setwd(temp_wd)
```

Then we use the `targets::tar_script()` function to write the `_targets.R` script with  the`target_Spec2Xtract()` factory:  

```{r}
## Write targets pipeline using the target_Spec2Xtract() factory
targets::tar_script(
  {
    require(Spec2Xtract)
    require(crew)

    ## The next option set the number of parallel workers to use
    tar_option_set(
      controller = crew_controller_local(workers = 2)
    )

    list(
      temp_tar <- target_Spec2Xtract(
        ## files should be the path(s) to the '.raw' files
        files = get_sample_rawfile(),
        ## cpd should be the compound table to extract (see the example)
        cpd = Spec2Xtract:::example_cpdlist_realdt,
        firstevent = TRUE,
        prec_ppm = 10,
        minscan = 3,
        rt_limit = 2,
        ppm = 10,
        save_dir = "./report"
      )
    )
  },
  ask = FALSE
)
```

Finally the pipeline can be run like this:  

```{r}
## Run the pipeline
targets::tar_make()
```

The results will be stored in the `./report` folder.

