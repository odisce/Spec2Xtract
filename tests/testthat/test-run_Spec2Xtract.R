test_that(
  "run_Spec2Xtract(): isopurity", {
    require(targets)
    raw_path <- normalizePath(
      file.path(
        testthat::test_path(),
        "testdata",
        "samplefile.raw"
      )
    )
    cpd_path <- normalizePath(
      file.path(
        testthat::test_path(),
        "testdata",
        "cpd.xlsx"
      )
    )
    targets::tar_dir({
      eval(substitute({
        run_Spec2Xtract(
          files_dir = dirname(raw_path),
          cpd_path = cpd_path,
          firstevent = TRUE,
          prec_ppm = 10,
          minscan = 3,
          rt_limit = 2,
          ppm = 10,
          ncore = 1,
          filter_irel = 0.5,
          filter_isopurity = 80,
          one_msp_file = TRUE,
          save_dir = "./report"
        )
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
        testthat::expect_true(any(grepl("CPD1", exported_files)))
        testthat::expect_true(any(grepl("File1", exported_files)))
        testthat::expect_true(any(grepl("SPEC", exported_files)))
        testthat::expect_true(any(grepl("spectra", exported_files)))
        testthat::expect_true(any(grepl("xics", exported_files)))
        testthat::expect_true(any(grepl("EventSummary", exported_files)))
        testthat::expect_true(any(grepl("Summary", exported_files)))
      },
      env = list(
        raw_path = raw_path,
        cpd_path = cpd_path
      )))
    })
  }
)
