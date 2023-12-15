## code to prepare `example_cpdlist` dataset goes here
example_cpdlist <- data.table::fread("./data-raw/cpd_list.txt")
usethis::use_data(example_cpdlist, overwrite = TRUE, internal = TRUE)
