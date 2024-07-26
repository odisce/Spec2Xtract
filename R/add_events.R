#' Get events info by parsing the scanType
#'
#' @param index_table The index table as returned by
#'                    [get_rawindex()]
#' @param firstevent Should all MSn events be checked individually
#'                   (`FALSE`) or only the first one (`TRUE`).
#'
#' @inheritParams check_rawfile
#'
#' @return
#' Return a `data.table` listing the different MS(n) events in
#' a **.raw** file. if `firstevent` is set to `TRUE`, the function
#' will check only the first event of each MSn events. If set to
#' `FALSE`, the function will check all events (it may take a lot
#' of time depending on the file).
#'
#' @export
#'
get_events_types <- function(
  index_table,
  rawpath,
  firstevent = TRUE
) {
  if (!"spec_energy" %in% names(index_table)) {
    index_table <- parse_index_dt(index_table = index_table)
  }
  event_table <- index_table[,
    .(scan_nb = .N),
    by = .(
      scanType,
      msLevel,
      spec_energy,
      spec_coltype,
      spec_polarity,
      spec_prec
    )
  ]

  if (firstevent == TRUE) {
    selected_events <- event_table[, utils::head(.SD, 1), by = msLevel]
    isowin_dt <- selected_events[, {
      scanType_sel <- scanType
      scan_to_check <- index_table[scanType == scanType_sel, ][1, scan]
      isowin <- rawrr::readSpectrum(
        rawfile = rawpath,
        scan = scan_to_check,
        mode = ""
      )[[1]]$`MS2 Isolation Width:`
      .(isolation_window = as.numeric(isowin))
    }, by = .(
      scanType,
      msLevel,
      spec_energy,
      spec_coltype,
      spec_polarity,
      spec_prec,
      scan_nb
    )]

    output <- merge(
      event_table,
      isowin_dt[, .(msLevel, isolation_window)],
      by = "msLevel"
    )
  } else {
    output <- event_table[, {
      scanType_sel <- scanType
      scan_to_check <- index_table[scanType == scanType_sel, ][1, scan]
      isowin <- rawrr::readSpectrum(
        rawfile = rawpath,
        scan = scan_to_check,
        mode = ""
      )[[1]]$`MS2 Isolation Width:`
      .(isolation_window = as.numeric(isowin))
    }, by = .(
      scanType,
      msLevel,
      spec_energy,
      spec_coltype,
      spec_polarity,
      spec_prec,
      scan_nb
    )]
  }
  return(output)
}
