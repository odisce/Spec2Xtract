test_that("Spec2Xtract::check_rawfile: .raw", {
  testthat::expect_true(Spec2Xtract::check_rawfile(rawrr::sampleFilePath()) == rawrr::sampleFilePath())
})
