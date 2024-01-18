#' Parse a spectrum object
#'
#' @param spectrum_list list of spectra
#' @param mz_range (optional) mass range to keep
#' @param mslevel Not sure DEV
#' @importFrom stringr str_squish
#' @import data.table magrittr
#'
#' @return
#' A list of spectra as `spectrum_list` but with 
#' formatted column and types.
#'
#' @export
#'
parse_spectrum_object <- function(
  spectrum_list,
  mz_range = NULL,
  mslevel = NULL
) {
  spectrums <- lapply(
    spectrum_list,
    function(spectrum_i) {
      mz_mode <- "NA"
      if (
        ("centroidStream" %in% names(spectrum_i)) &&
          isFALSE(spectrum_i$centroidStream)
      ) {
        if (
          ("HasCentroidStream" %in% names(spectrum_i)) &&
            grepl("True", spectrum_i$HasCentroidStream)
        ) {
          # centroid available
          mz_mode <- "centroid"
          mz_col <- "centroid.mZ"
          i_col <- "centroid.intensity"
        } else {
          # profile only
          mz_mode <- "profile"
          mz_col <- "mZ"
          i_col <- "intensity"
        }
      } else {
        warning("Unknown mz type")
        return(NULL)
      }
      spec_infos <- mapply(
        function(z, zz) {
          zz %>%
            {
              ifelse(
                . %in% names(spectrum_i),
                spectrum_i[[.]] %>% stringr::str_squish(),
                "NA"
              )
            }
        },
        c(
          "spec_scan",
          "spec_rt",
          "spec_iit",
          "spec_CID",
          "spec_HCDeV",
          "spec_HCD",
          "spec_isowin",
          "spec_etime",
          "spec_tic",
          "spec_mzrange",
          "spec_label",
          "spec_precmz"
        ),
        c(
          "scan",
          "rtinseconds",
          "Ion Injection Time (ms):",
          "API Source CID Energy:",
          "HCD Energy:",
          "HCD Energy eV:",
          "MS2 Isolation Width:",
          "Elapsed Scan Time (sec):",
          "TIC",
          "massRange",
          "scanType",
          "pepmass"
        ),
        SIMPLIFY = FALSE
      ) %>%
        as.list()

      spec_infos$spec_mslevel <- mslevel

      if (
        (
          (!is.null(mslevel) && mslevel > 1) ||
            is.null(mslevel)
        ) &&
          spec_infos$spec_label != "NA"
      ) {
        char_to_parse <- gsub(
          "^.*\\@([A-z]{1,6}[0-9]{1,4}\\.[0-9]{0,4}).\\[.*$",
          "\\1",
          x = spec_infos$spec_label
        )
        spec_infos$spec_energy <- gsub(
          "([A-z]{1,6})([0-9]{1,6}\\.[0-9]{0,6})",
          "\\2",
          char_to_parse
        )
        spec_infos$spec_colltype <- gsub(
          "([A-z]{1,6})([0-9]{1,6}\\.[0-9]{0,6})",
          "\\1",
          char_to_parse
        )
      } else {
        spec_infos$spec_energy <- as.character(NA)
        spec_infos$spec_colltype <- as.character(NA)
      }
      temp_i <- data.table::data.table(
        "mz" = spectrum_i[[mz_col]],
        "i" = spectrum_i[[i_col]]
      )
      if (!is.null(mz_range)) {
        temp_i <- temp_i[
          data.table::between(
            mz,
            mz_range[1],
            mz_range[2]
          )
        ]
      }
      output <- cbind(
        as.data.table(spec_infos),
        temp_i
      )
      return(output)
    }
  )
  return(spectrums)
}

#' Extract spectrum
#'
#' @param ms_scan_to_get scan index of the spectrum to extract
#' @inheritParams check_rawfile
#'
#' @return
#' Return a list with two levels:
#'   - `spectra_info_dt`: a data.table containing details on each
#'                        spectrum.
#'   - `spectra_db`: a data.table containing the spectrum (mz, i).
#'
#' @export
#'
#' @examples
#' get_spectrum_db(12, rawrr::sampleFilePath())
#'
get_spectrum_db <- function(ms_scan_to_get, rawpath) {
  scan_vi <- ms_scan_to_get[, unique(scan)]
  ms_spectrum <- rawrr::readSpectrum(
    rawpath,
    scan = scan_vi,
    mode = ""
  ) %>%
    {
      parse_spectrum_object(
        spectrum_list = .,
        mslevel = ms_scan_to_get[, unique(msLevel)]
      )
    }
  col_info <- ms_spectrum[[1]][,
    - c("mz", "i")
  ][,
    lapply(.SD, function(x) {length(unique(x))})
  ] %>%
    {
      unlist(.) %>%
        {
          which(. == 1)
        } %>%
        names()
    }
  spectrum_info <- ms_spectrum[[1]][, ..col_info] %>% unique()
  spectrum_info[, SpectrumIndex := ms_scan_to_get$SpectrumIndex]
  spectrum_info[, spec_polarity := ms_scan_to_get$spec_polarity]
  ms_spectrum_dt <- ms_spectrum[[1]][, .(mz, i)]
  return(
    list(
      "spectra_info_dt" = spectrum_info[],
      "spectra_db" = ms_spectrum_dt[]
    )
  )
}

#' Add spectrum to annobject
#'
#' @inheritParams add_events
#'
#' @import data.table
#'
#' @return
#' Fill the `spectra` slots of annobject.
#'
#' @export
#'
add_spectra <- function(annobject) {
  cpd_to_get <- sapply(annobject$cpd, function(x) {
    if (is.null(x$MSspectra$scan_info)) {
      return(FALSE)
    } else {
      nrow(x$MSspectra$scan_info) > 0
    }
  }) %>%
    which()
  for (i in cpd_to_get) {
    scan_dt <- annobject$cpd[[i]]$MSspectra$scan_info

    spectrum_db <- lapply(seq_len(nrow(scan_dt)), function(rowi) {
      scan_dt_i <- scan_dt[rowi, ]
      get_spectrum_db(
        ms_scan_to_get = scan_dt_i,
        rawpath = annobject$file$info$file_path[[scan_dt_i$FileIndex]]
      )
    })
    annobject$cpd[[i]]$MSspectra$spectra <- spectrum_db
  }
  return(annobject)
}
