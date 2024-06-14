test_that(
  "write_file: vector & overwrite", {
    tempa <- c("A", "B", "C")
    temp_file <- tempfile()
    res <- write_file(
      string = tempa,
      file_path = temp_file,
      overwrite = TRUE
    )
    expect_true(res)
    out <- readLines(temp_file)
    expect_true(all(out == tempa))

    res <- write_file(
      string = "D",
      file_path = temp_file,
      overwrite = FALSE
    )
    expect_true(res)
    out <- readLines(temp_file)
    expect_true(all(out == c(tempa, "D")))

    res <- write_file(
      string = "D",
      file_path = temp_file,
      overwrite = TRUE
    )
    expect_true(res)
    out <- readLines(temp_file)
    expect_true(all(out == "D"))
  }
)
