test_that("Spec2Xtract::add_events", {
  require(magrittr)
  require(data.table)
  temp_init <- Spec2Xtract::init_object(
    files = rawrr::sampleFilePath(),
    cpd = Spec2Xtract:::example_cpdlist_realdt
  )
  temp_event <- Spec2Xtract::add_events(
    temp_init,
    firstevent = TRUE
  )
  temp_events <- Spec2Xtract::add_cpd_events(temp_event, prec_ppm = 10)
  expect_true(data.table::is.data.table(temp_events$cpd[[1]]$MSEvents))
  expect_true(nrow(temp_events$cpd[[1]]$MSEvents) == 1)
  expect_true(
    all(
      c(
        "FileIndex",
        "EventIndex"
      ) %in% names(temp_events$cpd[[1]]$MSEvents)
    )
  )

  expect_error(Spec2Xtract::add_cpd_events(data.table()))
  temp_event_fail <- copy(temp_event)
  temp_event_fail$file$MSEvents <- NULL
  expect_error(Spec2Xtract::add_cpd_events(temp_event_fail))
  temp_event_fail <- copy(temp_event)
  temp_event_fail$cpd <- NULL
  expect_error(Spec2Xtract::add_cpd_events(temp_event_fail))
  temp_event_fail <- copy(temp_event)
  temp_event_fail$file$MSEvents[[1]] <- NULL
  temp_res <- Spec2Xtract::add_cpd_events(temp_event_fail)
  expect_true(is.data.table(temp_res$cpd[[1]]$MSEvents))
  expect_true(nrow(temp_res$cpd[[1]]$MSEvents) == 0)

  temp_add_events <- Spec2Xtract::cpd_add_events(
    cpd_dt = temp_events$cpd[[1]]$cpd_info,
    ms_events_dt = temp_events$file$MSEvents[[1]],
    prec_ppm = 10,
    indexonly = FALSE
  )
  expect_true(data.table::is.data.table(temp_add_events))
  expect_true(nrow(temp_add_events) == 1)
  expect_true(
    all(
      c(
        "spec_energy",
        "spec_coltype",
        "spec_polarity",
        "spec_prec"
      ) %in% names(temp_add_events)
    )
  )

  temp_add_events <- Spec2Xtract::cpd_add_events(
    cpd_dt = temp_events$cpd[[2]]$cpd_info,
    ms_events_dt = temp_events$file$MSEvents[[1]],
    prec_ppm = 10,
    indexonly = TRUE
  )
  expect_true(data.table::is.data.table(temp_add_events))
  expect_true(nrow(temp_add_events) == 1)
  expect_true(
    all(
      !c(
        "spec_energy",
        "spec_coltype",
        "spec_polarity",
        "spec_prec"
      ) %in% names(temp_add_events)
    )
  )


})
