test_that(
  "target_Spec2Xtract()", {
    require(targets)
    raw_path <- normalizePath(
      file.path(
        testthat::test_path(),
        "testdata",
        "samplefile.raw"
      )
    )
    targets::tar_dir({
      eval(
        substitute({
          dir.create("./report")
          targets::tar_script(
            {
              library(Spec2Xtract)
              library(targets)
              list(
                target_Spec2Xtract(
                  files = raw_path,
                  cpd = Spec2Xtract:::example_cpdlist_realdt,
                  firstevent = TRUE,
                  prec_ppm = 10,
                  minscan = 3,
                  rt_limit = 2,
                  ppm = 10,
                  save_dir = "./report"
                )
              )
            },
            ask = FALSE
          )
          ## Run pipeline
          targets::tar_make()
          ## Run tests
          targets::tar_read_raw("SUMMARY_dt") %>%
            {
              testthat::expect_true(
                is.data.table(.)
              )
            }

          targets::tar_read_raw("SUMMARYEVENT_dt") %>%
            {
              testthat::expect_true(
                is.data.table(.)
              )
            }

          targets::tar_read_raw("ISOPURITY") %>%
            {
              testthat::expect_true(
                is.list(.)
              )
            }

          exported_files <- list.files("./report", recursive = TRUE)
          expect_true(any(grepl("CPD1", exported_files)))
          expect_true(any(grepl("File1", exported_files)))
          expect_true(any(grepl("SPEC", exported_files)))
          expect_true(any(grepl("spectra", exported_files)))
          expect_true(any(grepl("xics", exported_files)))
          expect_true(any(grepl("EventSummary", exported_files)))
          expect_true(any(grepl("Summary", exported_files)))
        },
        env = list(raw_path = raw_path))
      )
    })
  }
)

testthat::test_that(
  "target_Spec2Xtract() with filter_irel and filter_isopurity", {
    require(targets)
    raw_path <- normalizePath(
      file.path(
        testthat::test_path(),
        "testdata",
        "samplefile.raw"
      )
    )
    filter_isopurity <- 80
    filter_irel <- 0.5
    targets::tar_dir({
      eval(
        substitute({
          dir.create("./report")
          targets::tar_script(
            {
              library(Spec2Xtract)
              library(targets)
              list(
                target_Spec2Xtract(
                  files = raw_path,
                  cpd = Spec2Xtract:::example_cpdlist_realdt,
                  firstevent = TRUE,
                  prec_ppm = 10,
                  minscan = 3,
                  rt_limit = 2,
                  ppm = 10,
                  filter_irel = filter_irel,
                  filter_isopurity = filter_isopurity,
                  save_dir = "./report"
                )
              )
            },
            ask = FALSE
          )
          ## Run pipeline
          targets::tar_make()
          ## Run tests
          targets::tar_read_raw("SPECTRA_DB") %>%
            {
              .$spectra_db
            } %>%
            {
              sapply(., function(x) {
                max(x$irel) >= filter_irel
              })
            } %>%
            {
              testthat::expect_true(
                all(.)
              )
            }

          exported_files <- list.files("./report", recursive = TRUE)

          ## Check that all exported files are with an isopurity above threshold
          event_noto_export <- targets::tar_read_raw("SUMMARYEVENT_dt") %>% {
            .[isopurity < filter_isopurity, SpecID]
          } %>%
            {
              sapply(
                .,
                function(x) {
                  any(grepl(x, grep("(xlsx|msp)$", exported_files, value = TRUE)))
                }
              )
            }
          testthat::expect_true(all(!event_noto_export))
          event_to_export <- targets::tar_read_raw("SUMMARYEVENT_dt") %>% {
            .[isopurity >= filter_isopurity, SpecID]
          } %>%
            {
              sapply(
                .,
                function(x) {
                  any(grepl(x, grep("(xlsx|msp)$", exported_files, value = TRUE)))
                }
              )
            }
          testthat::expect_true(all(event_to_export))
        },
        env = list(
          raw_path = raw_path,
          filter_isopurity = filter_isopurity,
          filter_irel = filter_irel
        ))
      )
    })
  }
)
