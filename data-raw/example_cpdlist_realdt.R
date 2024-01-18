## code to prepare `example_cpdlist` dataset goes here
example_cpdlist_realdt <- data.table::fread("./data-raw/cpd_list_realdt.txt")
usethis::use_data(example_cpdlist_realdt, overwrite = TRUE, internal = TRUE)
