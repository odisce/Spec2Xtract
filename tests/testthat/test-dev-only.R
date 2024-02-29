testthat::skip("For development purpose only, do not run !")

test_that("TODO", {
  # nolint start
  renv::update("Spec2Annot", prompt = FALSE)
  renv::snapshot(type = "explicit", prompt = FALSE)
  setwd("/home/sd265344/Spec2Xtract/")
  rm(list = ls())
  devtools::document()
  # devtools::test()
  devtools::check()
  renv::install(".", prompt = FALSE)
  # nolint end
})
