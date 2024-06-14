#' Wrapper to detect and extract spectra
#'
#' @inheritParams init_object
#' @inheritParams add_events
#' @inheritParams add_cpd_events
#' @inheritParams add_xics
#' @inheritParams add_best_peaks
#' @inheritParams add_events_to_extract
#' @inheritParams add_spectra
#' @inheritParams add_mspurity
#' @inheritParams add_annot
#' @inheritParams export_tables
#'
#' @return
#' If save_dir is not null, export the spectra and summary there.
#' Return a list with the following levels:
#'   - `"cpd_info"`: data.table with information on the compounds from cpd.
#'   - `"XICs"`: level to store XICs informations
#'   - `"Peaks"`: level to store peaks informations
#'   - `"MSspectra"`: level to store mass spectrum
#'
#' @export
#'
xtract_spectra <- function(
  files,
  cpd,
  firstevent = TRUE,
  prec_ppm = 10,
  ms1 = TRUE,
  rt_limit = 0.15,
  ppm = 5,
  save_dir = NULL,
  debug = FALSE
) {
  annoobj_test <- init_object(
    files = files,
    cpd = cpd
  ) %>%
    add_events(
      .,
      firstevent = firstevent
    ) %>%
    add_cpd_events(
      .,
      prec_ppm = prec_ppm
    ) %>%
    add_xics(
      .,
      ms1 = ms1,
      debug = debug
    ) %>%
    add_best_peaks(
      .,
      rt_limit = rt_limit,
      debug = debug
    ) %>%
    add_events_to_extract(
      .,
      debug = debug
    ) %>%
    add_spectra(.) %>%
    add_mspurity(
      .,
      debug = debug
    ) %>%
    add_annot(
      .,
      ppm = ppm
    )

  if (!is.null(save_dir)) {
    export_tables(
      annobject = annoobj_test,
      save_dir = save_dir
    )
  }

  return(annoobj_test)
}
