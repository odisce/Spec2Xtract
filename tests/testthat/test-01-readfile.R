test_that("Spec2Xtract::check_rawfile: .raw", {
  testthat::expect_true(
    Spec2Xtract::check_rawfile(rawrr::sampleFilePath()) == rawrr::sampleFilePath()
  )
})

test_that("Spec2Xtract::check_rawfile: .mzML", {
  f <- system.file("sciex/20171016_POOL_POS_1_105-134.mzML", package = "msdata")
  testthat::expect_true(
    Spec2Xtract::check_rawfile(f) == f
  )
})

