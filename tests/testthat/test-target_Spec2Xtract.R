test_that(
  "target_Spec2Xtract()", {
    require(targets)
    targets::tar_dir({
      dir.create("./report")
      targets::tar_script(
        {
          library(Spec2Xtract)
          library(targets)
          list(
            target_Spec2Xtract(
              files = file.path(
                testthat::test_path(),
                "testdata",
                "samplefile.raw"
              ),
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
      targets::tar_read_raw("SUMMARY_dt") %>% {
        testthat::expect_true(
          is.data.table(.)
        )
      }

      targets::tar_read_raw("SUMMARYEVENT_dt") %>% {
        testthat::expect_true(
          is.data.table(.)
        )
      }

      targets::tar_read_raw("ISOPURITY") %>% {
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
    })
  }
)
