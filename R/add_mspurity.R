#' Title
#'
#' @param spectrum_ms spectrum in data.table format with `mz` and `i`
#' @param msn_info details of the spectrum
#' @param cpd_info details of the compound (`mz_neg`, `mz_pos`)
#' @param prec_ppm mass tolerance to search the precursor (in ppm)
#' @importFrom Spec2Annot mz_range
#'
#' @return
#' A numeric value representing the precursor purity in the MS1
#' window (in percentage).
#' @export
#'
getpurity_from_spectrum <- function(
  spectrum_ms,
  msn_info,
  cpd_info,
  prec_ppm = 10
) {
  iso_win <- as.numeric(msn_info$spec_isowin)
  event_prec <- as.numeric(msn_info$spec_precmz)
  iso_range <- event_prec + c(-iso_win / 2, +iso_win / 2)
  polarity <- msn_info$spec_polarity %>% as.numeric()
  cpd_mz <- ifelse(
    polarity == 0,
    cpd_info$mz_neg,
    cpd_info$mz_pos
  )
  prec_range <- Spec2Annot::mz_range(cpd_mz, prec_ppm)
  spectrum_ms_iso <- spectrum_ms[mz %between% iso_range]
  prec_dt <- spectrum_ms_iso[
    mz %between% prec_range
  ][
    which.min(abs(mz - cpd_mz)), .(mz, i)
  ]
  isowin_purity <- prec_dt$i / spectrum_ms_iso[, sum(i)] * 100
  return(isowin_purity)
}

#' Add MS(n) purity
#'
#' @inheritParams add_events
#' @inheritParams add_xics
#'
#' @return
#' Fill the `IsoPurity` slot of annobject
#'
#' @export
#'
add_mspurity <- function(annobject, debug = FALSE) {
  cpd_to_get <- sapply(annobject$cpd, function(x) {
    if (is.null(x$MSspectra$spectra)) {
      return(FALSE)
    } else {
      length(x$MSspectra$spectra) > 0
    }
  }) %>%
    which()

  for (i in cpd_to_get) {
    if (isTRUE(debug)) {
      message("cpd: ", i, " file: ", appendLF = FALSE)
    }
    ## For each file
    for (
      file_i in annobject$cpd[[
        i
      ]]$MSspectra$scan_info[, unique(FileIndex)]
    ) {
      if (isTRUE(debug)) {
        message(file_i, "-", appendLF = FALSE)
      }
      scan_info_i <- annobject$cpd[[i]]$MSspectra$scan_info[FileIndex == file_i]
      spec_index_dt <- scan_info_i[msLevel == 2, ]
      if (nrow(spec_index_dt) > 0) {
        for (spec_ind_row in seq_len(nrow(spec_index_dt))) {
          if (isTRUE(debug)) {
            message(
              " event: ",
              spec_index_dt[spec_ind_row, SpectrumIndex],
              appendLF = FALSE
            )
          }
          spec_index_i <- spec_index_dt[spec_ind_row, ]
          msms_index <- spec_index_i$SpectrumIndex
          msms_info <- annobject$cpd[[i]]$MSspectra$spectra[[msms_index]]
          msms_info <- msms_info$spectra_info_dt
          if (
            scan_info_i[
              msLevel == 1 &
                spec_polarity == spec_index_i$spec_polarity,
              .N
            ] == 0
          ) {
            next
          }
          ms_index <- scan_info_i[
            msLevel == 1 &
              spec_polarity == spec_index_i$spec_polarity,
          ]
          si <- ms_index$SpectrumIndex

          isopurity_val <- getpurity_from_spectrum(
            spectrum_ms = annobject$cpd[[i]]$MSspectra$spectra[[si]]$spectra_db,
            msn_info = msms_info,
            cpd_info = annobject$cpd[[i]]$cpd_info,
            prec_ppm = 10
          )

          annobject$cpd[[i]]$MSspectra$spectra[[
            msms_index
          ]]$spectra_info_dt[, IsoPurity := isopurity_val]
        }
      }
    }
    if (isTRUE(debug)) {
      message("", appendLF = TRUE)
    }
  }
  return(annobject)
}
