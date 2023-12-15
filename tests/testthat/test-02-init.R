test_that("Spec2Xtract::get_rawindex: ", {
  temp_index <- Spec2Xtract::get_rawindex(rawrr::sampleFilePath())
  expect_true(data.table::is.data.table(temp_index))
  testthat::expect_true(
     == rawrr::sampleFilePath()
  )
})

test_that("Spec2Xtract::check_rawfile: .mzML", {
  f <- system.file("sciex/20171016_POOL_POS_1_105-134.mzML", package = "msdata")
  testthat::expect_true(
    Spec2Xtract::check_rawfile(f) == f
  )
})

