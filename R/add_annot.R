#' Annotate spectra
#'
#' @param ppm mass deviation for elemental composition
#'            search in ppm
#' @inheritParams add_events
#' @import data.table magrittr
#' @importFrom Spec2Annot annotate_mz
#'
#' @return
#' Annotate spectra using the compound elemental composition
#' with `Spec2Annot` and return an annotated annobject.
#'
#' @export
#'
add_annot <- function(annobject, ppm = 5) {
  cpd_to_get <- sapply(annobject$cpd, function(x) {
    if (is.null(x$MSspectra$spectra)) {
      return(FALSE)
    } else {
      length(x$MSspectra$spectra) > 0
    }
  }) %>%
    which()

  for (i in cpd_to_get) {
    cpd_scaninfo <- annobject$cpd[[i]]$MSspectra$scan_info
    msn_scan_dt <- cpd_scaninfo[msLevel > 1]
    if (nrow(msn_scan_dt) > 0) {
      for (scan_i in seq_len(nrow(msn_scan_dt))) {
        msms_index <- msn_scan_dt[scan_i, SpectrumIndex]
        msms_spectrum <- annobject$cpd[[
          i
        ]]$MSspectra$spectra[[msms_index]]$spectra_db[[1]][, .(mz, i)]

        cpd_compo <- annobject$cpd[[i]]$cpd_info$elemcomposition
        spec_annot <- Spec2Annot::annotate_mz(
          input_spectrum = msms_spectrum,
          ppm = ppm,
          polarity = msn_scan_dt[scan_i, spec_polarity],
          compo = cpd_compo,
          use_golden_ratio = FALSE
        )

        if (is.data.table(spec_annot)) {
          if (msms_spectrum[, .N] > spec_annot[, .N]) {
            spec_annot <- rbindlist(
              list(
                spec_annot,
                msms_spectrum[!ion_iter %in% spec_annot[, ion_iter]]
              ),
              fill = TRUE
            )
          }
          ## Save results
          annobject$cpd[[i]]$MSspectra$spectra[[
            msms_index
          ]]$spectra_db <- spec_annot[order(-irel), ]
        }
      }
    }
  }
  return(annobject)
}
