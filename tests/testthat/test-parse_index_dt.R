test_that(
  "parse_index_dt", {
    temp <- get_rawindex(rawrr::sampleFilePath())
    out <- parse_index_dt(temp)
    expect_true(is.data.table(out))
    expect_true(nrow(out) == nrow(temp))
    expect_true(
      all(
        c(
          "spec_energy",
          "spec_coltype",
          "spec_polarity",
          "spec_prec"
        ) %in% names(out)
      )
    )
    expect_true(all(1:2 %in% out[, unique(msLevel)]))
  }
)
