test_that(
  "get_sample_rawfile: url", {
    testthat::skip_if_offline()
    res <- get_sample_rawfile(url = NULL)
    expect_true(tools::file_ext(res) == "raw")
    expect_true(file.exists(res))
  }
)

test_that(
  "get_sample_rawfile: error", {
    temp <- get_sample_rawfile(url = "lmkqfmsdlkf")
    expect_false(temp)
  }
)
