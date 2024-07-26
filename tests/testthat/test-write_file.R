test_that(
  "write_file", {
    temp_path <- file.path(tempdir(), "test.msp")
    string_to_write <- c("mkljlmkjl", "lmkjlkj", "mklj")
    temp <- write_file(
      string = string_to_write,
      file_path = temp_path,
      overwrite = TRUE
    )
    expect_true(temp)
    expect_true(all(string_to_write %in% readLines(temp_path)))
    expect_true(length(readLines(temp_path)) == (length(string_to_write)))

    temp <- write_file(
      string = string_to_write,
      file_path = temp_path,
      overwrite = FALSE
    )
    expect_true(length(readLines(temp_path)) == (length(string_to_write) * 2))
    
    tempB <- write_file(
      string = c("mkljlmkjl", "lmkjlkj", "mklj"),
      file_path = temp_path,
      overwrite = TRUE
    )
    expect_true(length(readLines(temp_path)) == (length(string_to_write)))
  }
)
