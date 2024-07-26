

test_that(
  "open_cpd_file: file doesn't exist", {
    expect_error(
      open_cpd_file("/I/DONT/EXIST.tsv")
    )
  }
)

test_that(
  "open_cpd_file: wrong format", {
    cpd_path <- tempfile("test", fileext = ".lll")
    expect_error(open_cpd_file(cpd_path))
    write("lkj", cpd_path)
    expect_error(open_cpd_file(cpd_path))
  }
)

test_that(
  "open_cpd_file: existing file", {
    cpd_path <- test_get_path_cpd()
    cpd_dt <- open_cpd_file(cpd_path)
    ## Create .tsv format
    cpd_tsv_path <- tempfile(pattern = "test", fileext = ".tsv")
    fwrite(cpd_dt, file = cpd_tsv_path)
    
    lapply(
      c(cpd_tsv_path, cpd_path),
      function(x) {
        cpd_dt <- open_cpd_file(x)
        expect_true(is.data.table(cpd_dt))
        expect_true(nrow(cpd_dt) == 2)
        expect_true(
          all(
            c(
              "compound",
              "elemcomposition",
              "rtsec",
              "inchikey"
            ) %in% names(cpd_dt)
          )
        )
      }
    )
  }
)

test_that(
  "fun_check_cpd", {
    ## no inchikey
    ##
    cpd_dt <- test_get_path_cpd() %>%
      open_cpd_file(.)
    cpd_warn <- copy(cpd_dt)
    cpd_warn[, inchikey := NULL]
    cpd_err <- copy(cpd_dt)
    cpd_err[, rtsec := NULL]
    mapply(
      function(x, y) {
        if (y == 1) {
          expect_warning(fun_check_cpd(x))
        }
        if (y == 2) {
          expect_error(fun_check_cpd(x))
          return(TRUE)
        }
        cpd_dt_after_check <- fun_check_cpd(x)
        expect_true(all(names(x) %in% names(cpd_dt_after_check)))
        expect_true("CpdIndex" %in% names(cpd_dt_after_check))
        expect_true(nrow(x) == nrow(x))
      },
      list(cpd_dt, cpd_warn, cpd_err),
      c(0,1,2)
    )
  }
)

test_that(
  "fun_check_cpd: error", {
    expect_error(fun_check_cpd(data.table("tzeljk" = seq_len(4))))
  }
)