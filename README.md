
# Spec2Xtract <img src="man/figures/logo.svg" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/odisce/Spec2Xtract/branch/main/graph/badge.svg)](https://app.codecov.io/gh/odisce/Spec2Xtract?branch=main)
[![R-CMD-check](https://github.com/odisce/Spec2Xtract/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/odisce/Spec2Xtract/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Spec2Xtract** is an [R](https://www.r-project.org/) package to extract MS2 spectra from 
[*.mzML*](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#mzML) files or directly from *.raw* files 
acquired on [Thermo](https://www.thermofisher.com) instrument by using [mzR](https://github.com/sneumann/mzR/) 
or [rawrr](https://github.com/fgcz/rawrr) respectively. The package is compatible with multi-energy acquisition, 
extracting one spectra for every energy used in the file. Spectra annotation can be performed by using 
[**Spec2Annot**](https://github.com/odisce/Spec2Annot). An implementation on 
[Galaxy](https://workflow4metabolomics.usegalaxy.fr/) of those functionalities is under progress.

## Status

  [x] Main steps  
  [x] targets pipeline as a factory  
  [x] Documentation  
  [x] Tests  
  [ ] Vignette  

## Installation

You can install the development version of **Spec2Xtract** from [GitHub](https://github.com/) using your prefered method:

| method | command |
|--|--|
| [remotes](https://cran.r-project.org/web/packages/remotes/index.html) | `remotes::install_github("odisce/Spec2Xtract")` |
| [devtools](https://cran.r-project.org/web/packages/devtools/) | `devtools::install_github("odisce/Spec2Xtract")` |
| [renv](https://cran.r-project.org/web/packages/renv/) | `renv::install("github::odisce/Spec2Xtract")` |
| [pak](https://cran.r-project.org/web/packages/pak/) | `pak::pkg_install("odisce/Spec2Xtract")` |

**Spec2Xtract** use [**rawrr**](https://github.com/fgcz/rawrr) ro read *.raw* files directly. To be used it needs to 
be installed with the following commands:

```{r}
rawrr::installRawFileReaderDLLs()
rawrr::installRawrrExe()
```

## Inputs

To extract the MSn library, **Spec2Xtract** needs:  
  - a folder path (ex: */path/to/raw/directory/*) with the *.raw* or *.mzML* file(s). Any MSn acquisition is authorized in FIA mode or combined with LC or GC.  
  - a table with the neutral elemental composition (**elemcomposition**) of the molecules to extract (**compound**) and optionnaly the retention time (**rtsec**) (accepted formats *.txt*, *.csv*, *.tsv*, *.xlsx*). Any other column are optional and will be added in the final *.msp* library.  
    |  compound  | elemcomposition | rtsec |          inchikey           |
    | :--------: | :-------------: | :---: | :-------------------------: |
    |  Alanine   |     C3H7NO2     | 3.15  | QNAYBMKLOCPYGJ-REOHCLBHSA-N |
    | CustomName |    C4H21S1P3    | 1.52  |                             |


## Outputs

**Spec2Xtract** will extract every MSn events in range of a chromatographic peaks belonging to each of the molecules in the input table. It will automatically calculate the correct ions masses from the elemental composition with **Spec2Annot** and ouptut the results in a new folder *./report* with the following structure:  

```
└─ R-Project
  ├─ report
  │   ├─ figures
  |   |   ├─ spectras       # Contains all extracted spectra in image format
  |   |   |   ├─ CPD1_File1_SPEC1.png
  |   |   |   └─ ...
  |   |   └─ xics           # Contains all XICs with the detected peak windows
  |   |   |   ├─ CPD1.png
  |   |   |   └─ ...
  │   ├─ spectra
  |   |   ├─ msp            # Contains the library in .msp format
  |   |   |   └─ library.msp
  |   |   └─ xlsx           # Contains all extracted spectra in .xlsx format
  |   |   |   ├─ CPD1_File1_SPEC1.xlsx
  |   |   |   └─ ...
  │   ├─ EventSummary.xlsx  # Summary table of each MSn event extracted
  │   └─ Summary.xlsx       # Summary table of each initial compounds
  ├─ _targets               # targets cache
  └─ _targets.R             # targets pipeline
```

## Examples

### Using one wrapper

**Spec2Xtract** provide a wrapper which will initialize and run all the analysis with one function call:  
  - create a new project
  - open an R session inside
  - install **Spec2Xtract** and **rawrr** ([see Installation](#installation))
  - run the following command:  
      ```{r}
      Spec2Xtract::run_Spec2Xtract(
        files_dir = "/path/to/raw/directory/",
        cpd_path = "/path/to/cpd_info.xlsx",
        firstevent = TRUE,
        prec_ppm = 5,
        minscan = 3,
        rt_limit = 0.2,
        ppm = 6,
        save_dir = "./report",
        filter_irel = 0.01,
        filter_isopurity = 10,
        ncore = 1
      )
      ```
  - To get help: `?Spec2Xtract::run_Spec2Xtract`

### Using the targets factory implementation

For [targets](https://github.com/ropensci/targets) users, **Spec2Xtract** provide a [target factory](https://wlandau.github.io/targetopia/contributing.html) to run a full analysis pipeline
from a list of compound(s) and a list of path to *.raw* file(s) as shown below. [targets](https://github.com/ropensci/targets)
is a great package to manage analytical workflow and this implementation leverage all the functionalities it provides
(caching, vizualisation, CPUs parallelization, HPC deployment, etc...).

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


## to remove
xtract_spectra|init_object|add_events|add_cpd_events|add_xics|add_best_peaks|add_events_to_extract|add_spectra|add_mspurity|add_annot|export_tables|cpd_add_events