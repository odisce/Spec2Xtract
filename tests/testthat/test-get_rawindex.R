test_that(
  "get_rawindex", {
    temp <- get_rawindex(rawrr::sampleFilePath())
    expect_true(is.data.table(temp))
    expect_true(
      all(
        c(
          "scan",
          "scanType",
          "StartTime",
          "precursorMass",
          "MSOrder",
          "charge",
          "masterScan",
          "dependencyType",
          "monoisotopicMz"
        ) %in% names(temp)
      )
    )
    expect_true(nrow(temp) > 0)
  }
)