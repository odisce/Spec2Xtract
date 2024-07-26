test_that(
  "get_peak_range", {
    require(magrittr)
    temp_xic <- test_simulate_xic()
    get_peak_range(
      xic = temp_xic,
      direction = "L",
      peak_sel_center = temp_xic[which.max(i), scan]
    ) %>% {
      expect_true(. == 1)
    }
    get_peak_range(
      xic = temp_xic,
      direction = "R",
      peak_sel_center = temp_xic[which.max(i), scan]
    ) %>% {
      expect_true(. == temp_xic[, max(scan)])
    }
    get_peak_range(
      xic = temp_xic,
      direction = "N",
      peak_sel_center = temp_xic[which.max(i), scan]
    ) %>%
      expect_error()
  }
)

