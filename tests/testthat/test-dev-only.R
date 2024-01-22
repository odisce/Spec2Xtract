testthat::skip("For development purpose only, do not run !")

test_that("TODO", {
  ## Implementing:
  # init_object:           OK
  # add_events:            OK
  # add_cpd_events:        OK
  # add_xics:              OK
  # add_best_peaks:        OK to test
  # add_events_to_extract: OK to test
  # add_spectra:           OK to test
  # add_mspurity:          OK to test
  # add_annot:             OK to test
  # export_tables:         OK to test

  # renv::update("Spec2Annot", prompt = FALSE)
  # renv::snapshot(type = "explicit", prompt = FALSE)

  # devtools::document()
  # devtools::test()
  # renv::install(".", prompt = FALSE)

  # data_in <- exampleMS[10000:12000] ; SNR.Th <- 3
  # peakInfo <- peakDetectionCWT(data_in, SNR.Th = SNR.Th)

  # detect_centwave_peaks(data_in)

  # require(ggplot2)
  # temp_dt <- data.table(
  #   "scan" = seq_len(length(data_in)),
  #   "i" = data_in
  # )
  # ggplot(temp_dt, aes(scan, i)) +
  #   geom_line(data = peak_scale_trace_dt, linetype = 2, alpha = 0.75) +
  #   geom_line() +
  #   geom_line(
  #     data = temp_dt[scan %between% peak_info_dt[, c(scmin, scmax)]],
  #     color = "red"
  #   ) +
  #   geom_vline(
  #     xintercept = peak_info_dt$scapex,
  #     linetype = 2,
  #     color = "red"
  #   ) +
  #   geom_vline(
  #     xintercept = peak_info_dt[, c(scmin, scmax)], 
  #     linetype = 1,
  #     color = "red"
  #   ) +
  #   theme_bw()

  # majorPeakInfo <- peakInfo$majorPeakInfo
  # peakIndex <- majorPeakInfo$peakIndex
  # plotPeak(
  #   data_in,
  #   peakIndex,
  #   main = paste(
  #     "Identified peaks with SNR >",
  #     SNR.Th
  #   ),
  #   range = c(250, 700)
  # )
})


