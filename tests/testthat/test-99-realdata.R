file_path <- tryCatch(
  {
    url_sample <- "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS20/391_3-Hydroxy-4-methoxycinnamic_acid_NEG.RAW"
    temp_dir <- tempdir()
    temp_file <- file.path(temp_dir, "test.raw")
    download.file(url = url_sample, destfile = temp_file)
    temp_file
  },
  error = function(e) {
    print(e)
    return(FALSE)
  }
)

if (isFALSE(file_path) || !file.exists(file_path)) {
  testthat::skip(
    paste0(
      "Skip test: can't access sample file on:\n",
      "https://ftp.ebi.ac.uk/pub/databases/metabolights"
    )
  )
}

testthat::test_that("Real data", {
  require(data.table)
  require(magrittr)
  example_cpdlist_realdt <- Spec2Xtract:::example_cpdlist_realdt
  annoobj_test <- Spec2Xtract::init_object(
    files = file_path,
    cpd = example_cpdlist_realdt
  )

  annoobj_test <- Spec2Xtract::add_events(
    annobject = annoobj_test,
    firstevent = TRUE
  )

  annoobj_test <- Spec2Xtract::add_cpd_events(
    annobject = annoobj_test,
    prec_ppm = 3
  )

  annoobj_test <- Spec2Xtract::add_xics(
    annobject = annoobj_test,
    ms1 = TRUE,
    debug = TRUE
  )

  annoobj_test <- Spec2Xtract::add_best_peaks(
    annobject = annoobj_test,
    debug = TRUE
  )

  annoobj_test <- Spec2Xtract::add_events_to_extract(
    annobject = annoobj_test,
    debug = TRUE
  )

  annoobj_test <- Spec2Xtract::add_spectra(
    annobject = annoobj_test
  )

  annoobj_test <- Spec2Xtract::add_mspurity(
    annobject = annoobj_test,
    debug = TRUE
  )

   annoobj_test <- Spec2Xtract::add_annot(
        annobject = annoobj_test,
        ppm = 5
   )

  Spec2Xtract::export_tables(
    annobject = annoobj_test,
    save_dir = temp_dir
  )
  
})
