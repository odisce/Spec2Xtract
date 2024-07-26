test_that(
  "cpd_add_ionsmass", {
    require(magrittr)
    temp <- test_get_path_cpd() %>%
      open_cpd_file() %>%
      fun_check_cpd() %>%
      cpd_add_ionsmass()
    expect_true(is.data.table(temp))
    expect_true(
      all(
        c(
          "mz_pos",
          "mz_neg",
          "mz_neutral"
        ) %in% names(temp)
      )
    )
    expect_true(
      all(
        temp[, c(
          abs(mz_neutral - mz_pos),
          abs(mz_neutral - mz_neg)
        )] %between% c(1.007275, 1.007277)
      )
    )
  }
)
