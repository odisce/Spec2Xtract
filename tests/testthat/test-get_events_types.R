test_that(
  "get_events_types: firstevent = TRUE", {
    raw_path <- rawrr::sampleFilePath()
    raw_index <- get_rawindex(raw_path)
    raw_index <- raw_index[seq_len(5), ]
    output <- get_events_types(
      index_table = raw_index,
      rawpath = raw_path,
      firstevent = TRUE
    )
    expect_true(is.data.table(output))
    expect_true("isolation_window" %in% names(output))
    expect_true(is.numeric(output$isolation_window))
  }
)

test_that(
  "get_events_types: firstevent = FALSE", {
    raw_path <- rawrr::sampleFilePath()
    raw_index <- get_rawindex(raw_path)
    raw_index <- raw_index[seq_len(5), ]
    output <- get_events_types(
      index_table = raw_index,
      rawpath = raw_path,
      firstevent = FALSE
    )
    expect_true(is.data.table(output))
    expect_true("isolation_window" %in% names(output))
    expect_true(is.numeric(output$isolation_window))
  }
)
