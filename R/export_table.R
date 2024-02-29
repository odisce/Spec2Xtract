#' Title
#'
#' @inheritParams add_events
#'
#' @export
get_summary <- function(annobject) {
  output <- data.table()
  for (cpd_ind in seq_len(length(annobject$cpd))) {
    cpd_infos <- annobject$cpd[[cpd_ind]]$cpd_info
    cpd_infos[, FileNb := nrow(annobject$file$info)]

    # Peak Detection
    if (
      is.null(annobject$cpd[[cpd_ind]]$Peaks) ||
        nrow(annobject$cpd[[cpd_ind]]$Peaks) == 0
    ) {
      peakdetected <- 0
    } else {
      peakdetected <- nrow(annobject$cpd[[cpd_ind]]$Peaks)
    }

    cpd_infos[, PeakDetected := peakdetected]

    # MS2 Events
    if (
      is.null(
        annobject$cpd[[cpd_ind]]$MSspectra$scan_info
      ) ||
        nrow(
          annobject$cpd[[cpd_ind]]$MSspectra$scan_info[msLevel > 1]
        ) == 0
    ) {
      ms2events <- 0
    } else {
      ms2events <- nrow(
        annobject$cpd[[cpd_ind]]$MSspectra$scan_info[msLevel > 1]
      )
    }
    cpd_infos[, MS2Events := ms2events]

    ## Iso purity mean
    isopurity_mean <- sapply(
      annobject$cpd[[cpd_ind]]$MSspectra$spectra,
      function(x) {
        if (is.null(x)) {
          return(NULL)
        }
        if ("IsoPurity" %in% names(x$spectra_info_dt)) {
          return(x$spectra_info_dt$IsoPurity)
        } else {
          return(NULL)
        }
      }
    ) %>%
      unlist() %>%
      mean(., na.rm = TRUE)

    cpd_infos[, IsoPurity_mean := isopurity_mean]

    ## Spectra annotated
    if (!is.null(annobject$cpd[[cpd_ind]]$MSspectra$spectra)) {
      annotated_spec <- sapply(
        annobject$cpd[[cpd_ind]]$MSspectra$spectra,
        function(x) {
          if ("formula" %in% names(x$spectra_db)) {
            return(x$spectra_db[!is.na(formula), .N] > 0)
          } else {
            return(FALSE)
          }
        }
      ) %>%
        which() %>%
        length()
    } else {
      annotated_spec <- 0
    }
    cpd_infos[, annotated_spec := annotated_spec]
    cpd_infos[, cpd_index := cpd_ind]
    output <- rbind(
      cpd_infos,
      output
    )
  }
  setkey(output, "cpd_index")
  setcolorder(output)
  return(output[])
}

#' Title
#'
#' @inheritParams add_events
#'
#' @export
get_event_summary <- function(annobject) {
  output <- data.table()
  for (cpdi in seq_len(length(annobject$cpd))) {
    if (is.null(annobject$cpd[[cpdi]]$MSspectra)) {
      next
    }
    if (
      is.null(annobject$cpd[[cpdi]]$MSspectra$scan_info) ||
        annobject$cpd[[cpdi]]$MSspectra$scan_info[msLevel > 1, .N] == 0
    ) {
      next
    }
    spec_to_get <- annobject$cpd[[cpdi]]$MSspectra$scan_info[msLevel > 1, ]
    for (spec_rowi in seq_len(nrow(spec_to_get))) {
      spec_scans <- spec_to_get[spec_rowi, ]
      speci <- spec_scans[, SpectrumIndex]
      out_i <- annobject$cpd[[cpdi]]$MSspectra$spectra[[speci]]$spectra_info_dt
      out_i <- cbind(annobject$cpd[[cpdi]]$cpd_info, out_i[, -c("CpdIndex", "rtsec", "rtmin")])
      out_i[
        ,
        FileName := annobject$file$info[spec_scans[, FileIndex], file_name]
      ]
      out_i[, cpd_index := cpdi]
      out_i[
        ,
        spectrum_name := paste0(
          "CPD",
          cpdi,
          "_File",
          spec_scans$FileIndex,
          "_Spectrum",
          spec_scans$SpectrumIndex
        )
      ]

      output <- rbindlist(
        list(
          output,
          out_i[]
        ),
        fill = TRUE
      )
    }
  }

  if (nrow(output) <= 0) {
    return(FALSE)
  } else {
    setkey(output, "cpd_index", "spectrum_name")
    setcolorder(output)
    return(output[])
  }
}

#' Write spectrum to table
#'
#' @param save_dir Directory to save each spectrum
#' @inheritParams add_events
#'
#' @export
#'
export_spectrum_table <- function(annobject, save_dir) {
  for (cpdi in seq_len(length(annobject$cpd))) {
    if (!is.null(annobject$cpd[[cpdi]]$MSspectra)) {
      if (!is.null(annobject$cpd[[cpdi]]$MSspectra$spectra)) {
        if (
          !is.null(
            annobject$cpd[[cpdi]]$MSspectra$scan_info
          ) &&
            nrow(annobject$cpd[[cpdi]]$MSspectra$scan_info) > 0
        ) {
          if (
            annobject$cpd[[cpdi]]$MSspectra$scan_info[msLevel > 1, .N] > 0
          ) {
            event_to_get <- annobject$cpd[[
              cpdi
            ]]$MSspectra$scan_info[msLevel > 1, ]
            for (speci in seq_len(nrow(event_to_get))) {
              save_file <- event_to_get[
                speci,
                paste0(
                  "CPD",
                  cpdi,
                  "_File",
                  FileIndex,
                  "_Spectrum",
                  SpectrumIndex,
                  ".xlsx"
                )
              ]
              save_path <- file.path(save_dir, save_file)
              openxlsx::write.xlsx(
                x = annobject$cpd[[cpdi]]$MSspectra$spectra[[
                  event_to_get[speci, SpectrumIndex]
                ]]$spectra_db[order(-mz)],
                file = save_path
              )
            }
          }
        }
      }
    }
  }
}

#' Export spectra as tables
#'
#' @inheritParams add_events
#' @param save_dir Directory to save tables
#'
#' @import data.table openxlsx magrittr
#' @export
export_tables <- function(annobject, save_dir) {
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  save_l <- list(
    "summary_table" = file.path(save_dir, "Summary.xlsx"),
    "Events_table" = file.path(save_dir, "EventSummary.xlsx"),
    "Spectra_dir" = file.path(save_dir, "spectra")
  )
  save_l$Spectra_dir %>% {
    if (!dir.exists(.)) {
      dir.create(., recursive = TRUE)
    }
  }

  summary_table <- get_summary(annobject)
  openxlsx::write.xlsx(summary_table, file = save_l$summary_table)
  summary_events <- get_event_summary(annobject)

  if (!isFALSE(summary_events)) {
    openxlsx::write.xlsx(summary_events, file = save_l$Events_table)
  }

  # Get spectrums
  export_spectrum_table(annobject, save_dir = save_l$Spectra_dir)
  if (all(
    sapply(
      save_l,
      file.exists
    )
  )) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}