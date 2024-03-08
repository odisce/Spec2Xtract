test_that(
  "target_Spec2Xtract()", {
    require(targets)
    targets::tar_dir({
      path_data <- "./data"
      cpd_path <- file.path(path_data, "cpd.xlsx")
      dir.create(path_data)
      raw_path <- file.path(
        testthat::test_path(),
        "testdata",
        "samplefile.raw"
      )
      openxlsx::write.xlsx(
        Spec2Xtract:::example_cpdlist_realdt,
        cpd_path
      )
      file.copy(raw_path, file.path(path_data, "file.raw"))
      
      Spec2Xtract::run_Spec2Xtract(
        files_dir = path_data,
        cpd_path = cpd_path,
        firstevent = TRUE,
        prec_ppm = 10,
        minscan = 3,
        rt_limit = 2,
        ppm = 10,
        ncore = 1,
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
    })
  }
)
