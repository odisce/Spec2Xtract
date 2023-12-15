test_that("Spec2Xtract::get_rawindex", {
  temp_index <- Spec2Xtract::get_rawindex(rawrr::sampleFilePath())
  expect_true(data.table::is.data.table(temp_index))
  expect_true(all(dim(temp_index) == c(574,9)))
})

test_that("Spec2Xtract::parse_index_dt", {
  temp_index_parsed <- Spec2Xtract::parse_index_dt(index_table = Spec2Xtract::get_rawindex(rawrr::sampleFilePath()))
  expect_true(is.data.table(temp_index_parsed))
  expect_true(all(dim(temp_index_parsed) == c(574, 14)))
  expect_true(all(c("spec_energy", "spec_coltype", "spec_polarity", "spec_prec", "msLevel") %in% names(temp_index_parsed)))
  expect_true(temp_index_parsed[msLevel == 2 & is.na(spec_prec), .N] == 0)
  expect_true(temp_index_parsed[msLevel == 1 & is.na(spec_prec), .N] == temp_index_parsed[msLevel == 1, .N])
  expect_true(temp_index_parsed[msLevel == 2, unique(spec_coltype)] == "hcd")
})

# test_that("Spec2Xtract::fun_check_cpd", {
#   expect_true(TRUE)
#   Spec2Xtract::fun_check_cpd()
# })

# test_that('Spec2Xtract::cpd_add_ionsmass', {
#   expect_true(TRUE)
# })

# test_that('Spec2Xtract::init_object', {
#   expect_true(TRUE)
# })
