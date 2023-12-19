#' Get MS2 events for each compounds
#'
#' @param cpd_dt Compound table as returned by `Spec2Xtract::parse_index_dt()`
#' @param ms_events_dt data.table with MS events found in
#'                     `annobject$file$MSEvents`
#' @param prec_ppm Mass tolerance to find the precursor (in ppm)
#' @param indexonly Logical to return MS events info (`FALSE`) or ony
#'                  the index (`TRUE`)
#' @importFrom Spec2Annot mz_range
#' @return
#' A data.table containing the MS(n) events related to the compound.
#' If `indexonly = TRUE` only the MS(n) events indexes are returned,
#' else all the available informations are returned.
#' @export
#'
cpd_add_events <- function(
  cpd_dt,
  ms_events_dt,
  prec_ppm = 10,
  indexonly = FALSE
) {
  output <- cpd_dt[, {
    lapply(c(0, 1), function(x) {
      if (x == 0) {
        mz_range <- Spec2Annot::mz_range(mz_neg, prec_ppm)
      } else {
        mz_range <- Spec2Annot::mz_range(mz_pos, prec_ppm)
      }

      if (isTRUE(indexonly)) {
        ms_events_dt_iter <- rbind(
          ms_events_dt[spec_polarity == x & msLevel == 1, .(MSEvent_index)],
          ms_events_dt[
            spec_polarity == x &
              msLevel > 1 &
              (
                (spec_prec - (0.7 * isolation_window)) <= mz_range[1] &
                  (spec_prec + (0.7 * isolation_window)) >= mz_range[2]
              ),
            .(MSEvent_index)
          ]
        )
      } else {
        ms_events_dt_iter <- rbind(
          ms_events_dt[spec_polarity == x & msLevel == 1, ],
          ms_events_dt[
            spec_polarity == x & msLevel > 1 &
              (
                (spec_prec - (0.7 * isolation_window)) <= mz_range[1] &
                  (spec_prec + (0.7 * isolation_window)) >= mz_range[2]
              ),
          ]
        )
      }
      ms_events_dt_iter
    }) %>%
      rbindlist()
  }, by = names(cpd_dt)]
  return(output)
}

#' Add MSn events for each compounds
#'
#' @inheritParams add_events
#' @inheritParams cpd_add_events
#'
#' @return
#' Fill the `MSEvents` compound slot of the annobject
#'
#' @export
#'
add_cpd_events <- function(annobject, prec_ppm = 10) {
  if (!"list" %in% class(annobject)) {
    stop("annobject not recognized")
  }
  if (!"MSEvents" %in% names(annobject$file)) {
    stop("Need to run add_events() before this function")
  }
  if (!"cpd" %in% names(annobject)) {
    stop("cpd not initialized in object")
  }

  for (i in seq_len(length(annobject$cpd))) {
    if (isTRUE(annobject$cpd[[i]]$cpd_info$cpdCheck)) {
      ## Add events for each files
      annobject$cpd[[i]]$MSEvents <- lapply(
        seq_len(length(annobject$file$MSEvents)),
        function(x) {
          if (is.null(annobject$file$MSEvents[[x]])) {
            return(NULL)
          } else {
            cpd_add_events(
              cpd_dt = annobject$cpd[[i]]$cpd_info,
              ms_events_dt = annobject$file$MSEvents[[x]],
              prec_ppm = prec_ppm,
              indexonly = TRUE
            ) %>% {
              .[, .(FileIndex = x, MSEvent_index)]
            }
          }
        }
      ) %>%
        rbindlist()
    }
  }
  return(annobject)
}
