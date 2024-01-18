test_that("Spec2Xtract::get_events_types", {
  f <- rawrr::sampleFilePath()
  temp_events <- Spec2Xtract::get_events_types(
    index_table = Spec2Xtract::parse_index_dt(Spec2Xtract::get_rawindex(f)),
    raw_path = f,
    firstevent = TRUE
  )
  expect_true(data.table::is.data.table(temp_events))
  expect_true(
    all(
      c(
        "msLevel",
        "scanType",
        "spec_energy",
        "spec_coltype",
        "spec_polarity",
        "spec_prec",
        "scan_nb",
        "isolation_window"
      ) %in% names(temp_events)
    )
  )
  expect_true(all(c(1, 2) %in% temp_events[, unique(msLevel)]))
  expect_true("hcd" %in% temp_events[, unique(spec_coltype)])
  expect_true(all(!temp_events[, duplicated(spec_prec)]))

  temp_events_unique <- Spec2Xtract::get_events_types(
    index_table = Spec2Xtract::parse_index_dt(Spec2Xtract::get_rawindex(f)),
    raw_path = f,
    firstevent = FALSE
  )
  expect_true(all(c(1, 2) %in% temp_events_unique[, unique(msLevel)]))
  expect_true("hcd" %in% temp_events_unique[, unique(spec_coltype)])
  expect_true(all(!temp_events_unique[, duplicated(spec_prec)]))
})

test_that("Spec2Xtract::add_events", {
  temp_init <- Spec2Xtract::init_object(
    files = rawrr::sampleFilePath(),
    cpd = Spec2Xtract:::example_cpdlist_realdt
  )
  temp_event <- Spec2Xtract::add_events(temp_init, firstevent = TRUE)
  expect_true(
    data.table::is.data.table(
      temp_event$file$MSEvents[[1]]
    )
  )
  expect_true(
    nrow(
      temp_event$file$MSEvents[[1]]
    ) == temp_event$file$info$MSEvent_nb
  )
  expect_true(temp_event$file$info$Polarity == "1")
})
