#' Get events info by parsing the scanType
#'
#' @param index_table The index table as returned by
#'                    `Spec2Xtract::get_rawindex()`
#' @param raw_path Path to the raw file
#' @param firstevent Should all MSn events be checked individually
#'                        (`FALSE`) or only the first one (`TRUE`).
#'
#' @return
#' Return a `data.table` listing the different MS(n) events in
#' a **.raw** file. if `firstevent` is set to `TRUE`, the function
#' will check only the first event of each MSLevel. If set to
#' `FALSE`, the function will check all events (it may take a lot
#' of time depending on the file).
#'
#' @export
#'
get_events_types <- function(
  index_table,
  raw_path,
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
        rawfile = raw_path,
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
        rawfile = raw_path,
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

#' Add MSn events
#'
#' @param annobject Object build with `Spec2Extract::init_object()`
#' @inheritParams get_events_types
#'
#' @return
#' Implement the `MSEvents` and `info` slot of the object
#' @export
#'
add_events <- function(
  annobject,
  firstevent = TRUE
) {
  for (i in seq_len(nrow(annobject$file$info))) {
    if (annobject$file$info[i, FileCheck] == TRUE) {
      output <- get_events_types(
        index_table = annobject$file$index[[i]],
        raw_path = annobject$file$info[i, file_path],
        firstevent = firstevent
      )
      output[, EventIndex := seq_len(.N)]
      annobject$file$MSEvents[[i]] <- output[]
      ## Add info in fileinfo
      annobject$file$info[i,
        MSEvent_nb := output[, .N]
      ]
      annobject$file$info[i,
        Polarity := output[,
          paste0(unique(spec_polarity), collapse = "-") %>% as.character()
        ]
      ]
    } else {
      annobject$file$MSEvents[[i]] <- NULL
      annobject$file$info[i, MSEvent_nb := 0]
      annobject$file$info[i, Polarity := as.character(NA)]
    }
  }
  return(annobject)
}
