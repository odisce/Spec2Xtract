test_that("Spec2Xtract::check_rawfile: .raw", {
  testthat::expect_true(
    Spec2Xtract::check_rawfile(rawrr::sampleFilePath()) == rawrr::sampleFilePath()
  )

  f <- tempfile("test", fileext = ".raw")
  data.table::fwrite(x = data.table("test", "A", "B"), file = f)
  testthat::expect_error(
    Spec2Xtract::check_rawfile(f)
  )
})

test_that("Spec2Xtract::check_rawfile: .mzML", {
  f <- system.file("sciex/20171016_POOL_POS_1_105-134.mzML", package = "msdata")
  testthat::expect_true(
    Spec2Xtract::check_rawfile(f) == f
  )
  
  f <- tempfile("test", fileext = ".mzML")
  data.table::fwrite(x = data.table("test", "A", "B"), file = f)
  testthat::expect_error(
    Spec2Xtract::check_rawfile(f)
  )
})

test_that("Spec2Xtract::check_rawfile: .mzXML", {
  f <- system.file("threonine/threonine_i2_e35_pH_tree.mzXML", package = "msdata")
  testthat::expect_true(
    Spec2Xtract::check_rawfile(f) == f
  )
  
  f <- tempfile("test", fileext = ".mzXML")
  data.table::fwrite(x = data.table("test", "A", "B"), file = f)
  testthat::expect_error(
    Spec2Xtract::check_rawfile(f)
  )
})



test_that("Spec2Xtract::check_rawfile: Error file doesn't exist", {
  f <- "/sqdf/SQdfsq/DFsq/FDsqdF/sdqf.raw"
  testthat::expect_error(
    Spec2Xtract::check_rawfile(f)
  )
})

