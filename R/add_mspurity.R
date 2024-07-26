#' Title
#'
#' @param spectrum_ms spectrum in data.table format with `mz` and `i`
#' @param msn_info details of the spectrum with `spec_isowin`,
#'                 `spec_precmz` and `spec_polarity` columns
#' @param cpd_info details of the compound with `mz_neg` and `mz_pos`
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
  if (nrow(prec_dt) <= 0) {
    preci <- 0
  } else {
    preci <- prec_dt$i
  }
  isowin_purity <- preci / spectrum_ms_iso[, sum(i)] * 100
  return(isowin_purity)
}