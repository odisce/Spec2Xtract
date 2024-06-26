#' Get events
#'
#' @param pk_best data.table with details on the best peak
#' @inheritParams get_peak_scores_mult
#' @inheritParams get_events_types
#'
#' @return
#' a data.table with details on the MS(n) events to extract
#' for each compound.
#' @export
#'
get_events_scans <- function(cpd_events, pk_best, index_table) {
  if (!"CpdIndex" %in% names(cpd_events)) {
    if (
      ("CpdIndex" %in% names(pk_best)) &&
        (pk_best[, length(unique(CpdIndex))] == 1)
    ) {
      cpd_events[, CpdIndex := pk_best[, unique(CpdIndex)]]
    } else {
      cpd_events[, CpdIndex := 1]
      pk_best[, CpdIndex := 1]
    }
  }

  output <-  cpd_events[, {
    CpdIndex_i <- unique(CpdIndex)
    pk_info <- pk_best[CpdIndex == CpdIndex_i, ]
    sub_dt <- .SD[
      ,
      .(
        msLevel,
        spec_energy,
        spec_coltype,
        spec_polarity,
        scanType
      )
    ]
    raw_sub <- index_table[
      sub_dt,
      on = c(
        "msLevel",
        "spec_energy",
        "spec_coltype",
        "spec_polarity",
        "scanType"
      )
    ]
    ms2_scan <- raw_sub[
      data.table::between(
        StartTime,
        pk_info$rtmin,
        pk_info$rtmax
      ),
    ]
    ## Select closest scan for each type of events
    ms2_scan <- ms2_scan[,
      .SD[
        which.min(
          abs(StartTime - (pk_info$rt))
        ),
      ],
      by = .(
        msLevel,
        spec_energy,
        spec_coltype,
        spec_polarity
      )
    ]
    ms2_scan[, .(
      scan,
      msLevel,
      spec_energy,
      spec_coltype,
      spec_polarity,
      spec_prec,
      scanType,
      rtsec = as.numeric(StartTime) * 60,
      rtmin = as.numeric(StartTime),
      masterScan = as.integer(masterScan)
    )]
  }, by = .(CpdIndex)]
  return(output)
}

#' Add MS(n) events to extract
#'
#' @inheritParams add_events
#' @inheritParams add_xics
#'
#' @return
#' Fill the `MSspectra` slots of the annobject with
#' information on MS(n) event to extract in the next
#' step.
#'
#' @export
#'
add_events_to_extract <- function(annobject, debug = FALSE) {
  cpd_to_get <- sapply(annobject$cpd, function(x) {
    nrow(x$Peaks) > 0
  }) %>%
    which()
  for (i in cpd_to_get) {
    if (isTRUE(debug)) {
      message("cpd: ", i, " file: ", appendLF = FALSE)
    }
    output_dt <- data.table()
    for (filei in annobject$cpd[[i]]$Peaks$FileIndex) {
      if (isTRUE(debug)) {
        message(filei, "-", appendLF = FALSE)
      }
      event_to_extract_dt <- get_events_scans(
        cpd_events = annobject$file$MSEvents[[filei]][
          EventIndex %in%
            annobject$cpd[[i]]$MSEvents[
              FileIndex == filei, EventIndex
            ]
        ],
        pk_best = annobject$cpd[[i]]$Peaks[FileIndex == filei],
        index_table = annobject$file$index[[filei]]
      )
      if (nrow(event_to_extract_dt) > 0) {
        event_to_extract_dt[, FileIndex := filei]
        output_dt <- rbind(event_to_extract_dt, output_dt)
      } else {
        next
      }
    }
    if (isTRUE(debug)) {
      message("", appendLF = TRUE)
    }
    output_dt[, SpectrumIndex := seq_len(.N)]
    annobject$cpd[[i]]$MSspectra$scan_info <- output_dt[]
  }

  return(annobject)
}
