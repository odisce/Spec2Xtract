test_that(
  "get_peaks_xic", {
    xics <- test_get_trueXIC()
    output <- get_peaks_xic(
      xics,
      minscan = 2,
      majpeakonly = FALSE
    )
    expect_true(is.data.table(output))
    expect_true(nrow(output) == 1)
    expect_true(
      all(
        c(
          "peakID",
          "scmin",
          "scmax",
          "scpos",
          "scale",
          "into",
          "maxo",
          "sn",
          "rtmin",
          "rtmax",
          "rt",
          "cpd_index",
          "scan_nb",
          "zigzag_score"
        ) %in% names(output)
      )
    )
    expect_warning(
      get_peaks_xic(
        xics,
        minscan = 2,
        majpeakonly = TRUE
      )
    )
  }
)
