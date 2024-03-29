test_that("Spec2Xtract::get_rawindex", {
  temp_index <- Spec2Xtract::get_rawindex(rawrr::sampleFilePath())
  expect_true(data.table::is.data.table(temp_index))
  expect_true(all(dim(temp_index) == c(574, 9)))
})

test_that("Spec2Xtract::parse_index_dt", {
  temp_index_table <- Spec2Xtract::get_rawindex(rawrr::sampleFilePath())
  temp_index_parsed <- Spec2Xtract::parse_index_dt(
    index_table = temp_index_table
  )
  expect_true(is.data.table(temp_index_parsed))
  expect_true(all(dim(temp_index_parsed) == c(574, 14)))
  expect_true(
    all(
      c(
        "spec_energy",
        "spec_coltype",
        "spec_polarity",
        "spec_prec",
        "msLevel"
      ) %in% names(temp_index_parsed)
    )
  )
  expect_true(
    temp_index_parsed[msLevel == 2 & is.na(spec_prec), .N] == 0
  )
  expect_true(
    temp_index_parsed[
      msLevel == 1 & is.na(spec_prec), .N
    ] == temp_index_parsed[msLevel == 1, .N]
  )
  expect_true(
    temp_index_parsed[
      msLevel == 2, unique(spec_coltype)
    ] == "hcd"
  )
})

test_that("Spec2Xtract::fun_check_cpd", {
  temp_cpd <- Spec2Xtract::fun_check_cpd(Spec2Xtract:::example_cpdlist_realdt)
  expect_true(is.data.table(temp_cpd))
  expect_true(
    all(
      c(
        "compound",
        "elemcomposition",
        "rtsec",
        "inchikey",
        "CpdIndex",
        "rtmin"
      ) %in% names(temp_cpd)
    )
  )
  expect_true(temp_cpd[, all(is.na(inchikey))])
  expect_true(temp_cpd[, all(is.numeric(rtsec))])
  expect_true(temp_cpd[, all(is.numeric(rtmin))])
  expect_true(temp_cpd[, all(is.character(compound))])
  expect_true(temp_cpd[, all(is.character(elemcomposition))])
  expect_true(temp_cpd[, all(is.integer(CpdIndex))])
  ## Error
  expect_error(
    Spec2Xtract::fun_check_cpd(
      data.table("qldsjkf", "compound")
    )
  )
  expect_warning(
    Spec2Xtract::fun_check_cpd(
      data.table(
        "rtsec" = 10.5,
        "compound" = "A",
        "elemcomposition" = "C6H12"
      )
    )
  )
})

test_that("Spec2Xtract::cpd_add_ionsmass", {
  temp_cpd <- Spec2Xtract::fun_check_cpd(Spec2Xtract:::example_cpdlist_realdt)
  temp_cpd_supp <- Spec2Xtract::cpd_add_ionsmass(temp_cpd)
  expect_true(is.data.table(temp_cpd_supp))
  expect_true(all(dim(temp_cpd_supp) == dim(temp_cpd) + c(0, 3)))
  expect_true(
    all(
      c(
        "mz_neutral",
        "mz_pos",
        "mz_neg"
      ) %in% names(temp_cpd_supp)
    )
  )
  expect_true(
    all(
      temp_cpd_supp[,
        lapply(.SD, is.numeric),
        .SDcols = c(
          "mz_neutral",
          "mz_pos",
          "mz_neg"
        )
      ]
    )
  )
  expect_true(
    all(
      temp_cpd_supp[, mz_pos > mz_neutral & mz_neg < mz_neutral]
    )
  )
})

test_that("Spec2Xtract::init_object", {
  temp_init <- Spec2Xtract::init_object(
    files = rawrr::sampleFilePath(),
    cpd = Spec2Xtract:::example_cpdlist_realdt
  )
  expect_true(is.list(temp_init))
  expect_true(all(c("file", "cpd") %in% names(temp_init)))
  expect_true(is.list(temp_init$file))
  expect_true(all(c("info", "index") %in% names(temp_init$file)))
  expect_true(nrow(temp_init$file$info) == 1)
  expect_true(length(temp_init$file$index) == 1)
  expect_true(length(temp_init$cpd) == nrow(Spec2Xtract:::example_cpdlist_realdt))
  expect_true("cpd_info" %in% names(temp_init$cpd[[1]]))
  expect_true(is.data.table(temp_init$cpd[[1]]$cpd_info))
  expect_true(temp_init$cpd[[1]]$cpd_info$CpdIndex == 1)
  expect_true(temp_init$cpd[[2]]$cpd_info$CpdIndex == 2)
  expect_true("cpdCheck" %in% names(temp_init$cpd[[2]]$cpd_info))
  expect_true(temp_init$cpd[[2]]$cpd_info$cpdCheck == TRUE)
  ## Error file
  f <- tempfile("test", fileext = ".raw")
  data.table::fwrite(x = data.table("A", "C", "D"), file = f)
  temp_init <- Spec2Xtract::init_object(
    files = f,
    cpd = Spec2Xtract:::example_cpdlist_realdt
  )
  expect_true(temp_init$file$info$FileCheck == FALSE)
  ## Error cpd
  cpd_dt <- copy(Spec2Xtract:::example_cpdlist_realdt)
  cpd_dt[1, elemcomposition := "NA150-lkjSDF10"]
  temp_init <- Spec2Xtract::init_object(
    files = rawrr::sampleFilePath(),
    cpd = cpd_dt
  )
  expect_true(temp_init$cpd[[1]]$cpd_info$cpdCheck == FALSE)
})
