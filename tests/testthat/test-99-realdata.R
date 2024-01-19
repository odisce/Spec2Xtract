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
  save_dir <- tempdir()

  temp_res <- xtract_spectra(
    files = file_path,
    cpd = example_cpdlist_realdt,
    firstevent = TRUE,
    prec_ppm = 10,
    ms1 = TRUE,
    rt_limit = 1,
    ppm = 5,
    save_dir = save_dir,
    debug = FALSE
  )

  expect_true(is.list(temp_res))
  expect_true(all(c("file", "cpd") %in% names(temp_res)))
})
