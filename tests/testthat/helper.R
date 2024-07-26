test_get_path_raw <- function() {
  raw_path <- normalizePath(
    file.path(
      testthat::test_path(),
      "testdata",
      "samplefile.raw"
    )
  )
  return(raw_path)
}

test_get_path_cpd <- function() {
  cpd_path <- normalizePath(
    file.path(
      testthat::test_path(),
      "testdata",
      "cpd.xlsx"
    )
  )
  return(cpd_path)
}

test_simulate_xic <- function() {
  require(magrittr)
  temp_xic <- seq(1, 100, by = 10) %>% {
    c(., rev(.)[-1])
  } %>% {
    data.table("i" = .)
  }
  temp_xic[, scan := seq_len(.N)]
  temp_xic[, rt := seq_len(.N)]
  return(temp_xic[])
}

test_get_trueXIC <- function() {
  require(magrittr)
  require(data.table)
  fun_get_xic(
    rawpath = test_get_path_raw(),
    mz = test_get_path_cpd() %>%
      open_cpd_file() %>%
      fun_check_cpd() %>%
      cpd_add_ionsmass() %>% {
        .[1, mz_neg]
      },
    ppm = 5,
    filter = "Ms"
  )
}
