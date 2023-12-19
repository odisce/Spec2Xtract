test_that("Spec2Xtract::add_xics", {
  temp_init <- Spec2Xtract::init_object(
    files = rawrr::sampleFilePath(),
    cpd = Spec2Xtract:::example_cpdlist
  )
  temp_event <- Spec2Xtract::add_events(
    temp_init,
    firstevent = TRUE
  )
  temp_events <- Spec2Xtract::add_cpd_events(temp_event, prec_ppm = 10)

  temp <- Spec2Xtract::add_xics(
    annobject = temp_events,
    ms1 = TRUE,
    debug = TRUE
  )

  expect_true(data.table::is.data.table(temp$cpd[[1]]$XICs))

  expect_true(
    all(
      c(
        "ppm",
        "rt",
        "i",
        "type",
        "MSEvent_index",
        "FileIndex"
      ) %in%
        names(temp$cpd[[1]]$XICs)
    )
  )
})
