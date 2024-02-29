#' Calculate peak scores
#'
#' @param peaks_dt Peak list as returned by
#'                  `Spec2Xtract::get_peaks_from_xic()`
#' @param ref_rtmin reference rt in minute
#' @import data.table magrittr
#' @importFrom stats na.omit
#' @return
#' A data.table identical to `peaks_dt` with the following
#' supplementary columns:
#'   - `peak_in_range`: Number of peaks in range for this MS Event
#'   - `diff_from_ref`: RT shift from reference
#'   - `peak_score`: Final score
#' @export
#'
get_peak_scores <- function(
  peaks_dt,
  ref_rtmin = NULL
) {
  peaks_dt[, "iterPeaks" := seq_len(.N)]
  peaks_dt[
    rtmin < ref_rtmin & rtmax > ref_rtmin,
    "peak_in_range" := .N,
    by = filter
  ]
  peaks_dt[, "diff_from_ref" := abs(rt - ref_rtmin)]
  max_into <- peaks_dt[, max(into, na.rm = TRUE)]
  peaks_dt[
    ,
    "peak_score" := c(
      (1 / zigzag_score),
      ((into / max_into) * 15),
      (1 / diff_from_ref * 2),
      peak_in_range * 3,
      ifelse(scan_nb < 6, scan_nb, 10),
        ifelse(sn > 100, 10,
          ifelse(sn > 50, 5,
            ifelse(sn > 10, 10, 1)
          )
        )
    ) %>%
      na.omit() %>%
      sum(),
    by = .(iterPeaks)
  ]
  peaks_dt[, iterPeaks := NULL]
  return(peaks_dt[])
}

#' Calculate peak scores
#'
#' @param xic_peaks Peak list as returned by
#'                  `Spec2Xtract::get_peaks_from_xic()`
#' @param cpd_events Compound MS(n) event slot
#' @param CpdIndex (optional) Compound iteration ID to use
#' @param cpd_info Compound Info slot or reference rt in minute
#' @import data.table magrittr
#' @importFrom stats na.omit
#' @return
#' The `xic_peaks` data.table with a new columns:
#'   - `peak_score`: Overall scoring of the peaks
#'   - `peak_in_range`: Number of peaks in range
#'   - `peaks_mslevels`: MS(n) level of the peak
#'   - `diff_from_ref`: Time difference from reference
#' @export
#'
get_peak_scores_mult <- function(xic_peaks, cpd_events, CpdIndex = NULL, cpd_info = NULL) {
  if (!is.null(CpdIndex)) {
    if (CpdIndex %in% xic_peaks[, unique(CpdIndex)]) {
      CpdIndex_sel <- CpdIndex
    } else {
      stop("CpdIndex not found in xic_peaks$CpdIndex")
    }
  } else {
    if (!"CpdIndex" %in% names(xic_peaks)) {
      xic_peaks[, CpdIndex := 1]
    }
    CpdIndex_sel <- xic_peaks[, unique(CpdIndex)]
  }
  output <- lapply(CpdIndex_sel, function(x) {
    cpd_dt_i <- xic_peaks[CpdIndex == x]
    peak_scores_i <- get_peak_scores(cpd_dt_i, cpd_info[CpdIndex == x, rtmin])
    return(peak_scores_i[])
  }) %>%
    rbindlist()
  return(output)
}

#' Use MassSpecWavelet to detect peaks
#'
#' @param xic The XIC as a data.table with the columns:
#'            'scan' and 'i'.
#' @param direction Character to choose the right 'R' or left 'L'
#'                  side of the peak
#' @param peak_sel_center Scan corresponding to the peak apex.
#' @import data.table
#' @return
#' A numeric value corresponding to the scan on the right `R` or
#' left `L` side of the peak.
#'
#' @export
#'
get_peak_range <- function(xic, direction = c("R", "L"), peak_sel_center) {
  temp_xic <- copy(xic)
  if (direction == "R") {
    temp_xic[, i_n := i - shift(i, 1)]
    temp_side <- temp_xic[scan > peak_sel_center]
  } else if (direction == "L") {
    temp_xic[, i_n := i - shift(i, -1)]
    temp_side <- temp_xic[scan < peak_sel_center]
  } else {
    stop("direction arg not one of 'L' or 'R'")
  }

  if (nrow(temp_side[i_n >= 0, ]) > 0) {
    output <- switch(
      direction,
      "R" = {
        temp_side[i_n >= 0, ][, min(scan)-1]
      },
      "L" = {
        temp_side[i_n >= 0, ][, max(scan)+1]
      },
      stop("direction arg not one of 'L' or 'R'")
    )
    return(output)
  } else {
    output <- switch(
      direction,
      "R" = {
        temp_side[, max(scan)]
      },
      "L" = {
        temp_side[, min(scan)]
      },
      stop("direction arg not one of 'L' or 'R'")
    )
    return(output)
  }
}

#' Use MassSpecWavelet to detect peaks
#'
#' @param vecint Numeric vector containing the intensities
#'               sorted by the time dimension.
#' @param snrth Signal over Noise threshold (see
#'               `MassSpecWavelet::peakDetectionCWT()`)
#' @param majpeakonly Logical to extract only the majors peaks (`TRUE`)
#'                    or all peaks (`FALSE`) from MassSpecWavelet
#' @param limits Should peak boundaries be found on true data (`xic`) or
#'               on CentWave scales (`scales`: default)?
#' @param ... Other options passed to `MassSpecWavelet::peakDetectionCWT()`
#' @importFrom MassSpecWavelet peakDetectionCWT prepareWavelets
#' @return
#' A data.table with information on chromatographic peaks
#' with the following columns:
#'   - `peakID`: Peak name given by `MassSpecWavelet::peakDetectionCWT()`
#'   - `scmin`, `scmax`, `scpos`: left, right boundaries and apex of the peak
#'                                (index of the values in vecint)
#'   - `scale`: scale of the peak given by `MassSpecWavelet::peakDetectionCWT()`
#' @export
#'
detect_centwave_peaks <- function(
  vecint,
  snrth = 3,
  majpeakonly = TRUE,
  limits = c("xic", "scales")[2],
  ...
) {
  vecint <- c(0, 0, 0, vecint, 0, 0, 0)
  peak_info <- MassSpecWavelet::peakDetectionCWT(
    vecint,
    scales = MassSpecWavelet::prepareWavelets(
      scales = c(seq(1, 10, 1), seq(10, 30, 2), seq(32, 64, 4)),
      mslength = length(vecint),
      wavelet_xlimit = 3
    ),
    SNR.Th = snrth,
    tuneIn = FALSE,
    ...
  )
  if (isTRUE(majpeakonly)) {
    peak_to_get <- names(peak_info$majorPeakInfo$peakIndex)
  } else {
    peak_to_get <- names(peak_info$majorPeakInfo$allPeakIndex)
  }
  if (length(peak_to_get) <= 0) {
    warning(
      paste0(
        "detect_centwave_peaks(): ",
        "No peak found with the selected parameters"
      )
    )
    return(FALSE)
  } else {
    output <- lapply(
      peak_to_get,
      function(peak_sel_name) {
        peak_sel_apex <-
          peak_info$majorPeakInfo$allPeakIndex[[peak_sel_name]]
        peak_sel_center <-
          peak_info$majorPeakInfo$peakCenterIndex[[peak_sel_name]]
        peak_sel_scale <-
          peak_info$majorPeakInfo$peakScale[[peak_sel_name]]
        peak_sel_sn <-
          peak_info$majorPeakInfo$peakSNR[[peak_sel_name]]
        ## Find rt boundaries
        xcoefs_num <- as.numeric(colnames(peak_info$wCoefs))
        
        ## using scales or true xics
        if (limits == "xic") {
          peak_scale_trace <- vecint
        } else if (limits == "scales") {
          peak_scale_trace <- peak_info$wCoefs[,
            which.min(abs(xcoefs_num - peak_sel_scale))
          ]
        } else {
          stop("limits argument not recognized, should be one of 'xic' or 'scales'")
        }


        peak_scale_trace_dt <- data.table(
          "scan" = seq_len(length(peak_scale_trace)),
          "i" = peak_scale_trace
        )

        # peak_scale_trace_dt[, i_n := i - shift(i, 1)]
        # peak_scale_trace_dt[, i_p := i - shift(i, -1)]
        scmin <- get_peak_range(
          peak_scale_trace_dt,
          "L",
          peak_sel_center
        ) - 3
        scmax <- get_peak_range(
          peak_scale_trace_dt,
          "R",
          peak_sel_center
        ) - 3
        scmin <- ifelse(scmin < 1, 1, scmin)
        scmax <- ifelse(scmax > (length(vecint) - 6), length(vecint) - 6, scmax)
        scpos <- peak_sel_apex - 3
        if (scpos <= 0) {
          scpos <- 1
        } else if (scpos > (length(vecint) - 6)) {
          scpos <- (length(vecint) - 6)
        }

        peak_info_dt <- data.table(
          "peakID" = peak_sel_name %>% as.character(),
          "scmin" = scmin %>% as.integer(),
          "scmax" = scmax %>% as.integer(),
          "scpos" = scpos %>% as.integer(),
          "scale" = peak_sel_scale %>% as.integer(),
          "into" = sum(vecint[scmin:scmax], na.rm = TRUE) %>% as.numeric(),
          "maxo" = max(vecint[scmin:scmax], na.rm = TRUE) %>% as.numeric(),
          "sn" = peak_sel_sn %>% as.numeric()
        )
        return(peak_info_dt)
      }
    )
    data.table::rbindlist(output, fill = TRUE)
  }
}

#' Calculate the zig-zag score of a vector
#'
#' @param vecint Numeric vector containing the intensities
#'               sorted by the time dimension.
#'
#' @return The zig-zag index value
#' @source \url{https://github.com/GauravPandeyLab/MetaClean/blob/master/R/calculateZigZagIndex.R}
#' @export
#'
#' @examples
#' calculate_zigzag(rnorm(100, 1e5, sd = 1e3))
#' calculate_zigzag(rnorm(100, 1e5, sd = 1e3))
calculate_zigzag <- function(vecint) {
  if (length(vecint) > 4) {
    eic <- vecint
    end <- length(eic)
    epi <- max(eic) - (eic[1] + eic[2] + eic[end] + eic[end - 1]) / 4
    zig_zag_sum <- 0
    for (i in 2:(end - 1)) {
      local_zig_zag <- (2 * eic[i] - eic[i - 1] - eic[i + 1])^2
      zig_zag_sum <- zig_zag_sum + local_zig_zag
    }
    zig_zag_index <- zig_zag_sum / (epi^2 * end)
  } else {
    zig_zag_index <- NA
  }
  return(zig_zag_index)
}

#' Find peaks usint CentWave
#'
#' @param xic_mat a matrix or data.table with the following columns:
#'                  - `rt` with the retention time
#'                  - `i` with the intensities
#' @param minscan numeric value for the minimum number of scan
#'                inside a peak to keep it
#' @inheritParams detect_centwave_peaks
#'
#' @export
#' @import data.table
#'
get_peaks_xic <- function(xic_mat, minscan = 3, majpeakonly = FALSE) {
  if (data.table::is.data.table(xic_mat)) {
    if (any(!c("rt", "i") %in% names(xic_mat))) {
      stop("if xic_mat is a data.table, we need 'rt' and 'i' columns")
    }
  } else if (is.matrix(xic_mat)) {
    if (any(!c("i") %in% colnames(xic_mat))) {
      stop("if xic_mat is a matrix, we need a i column")
    }
  } else {
    stop(
      paste0(
        "Input not recognized, must be a data.table (with ",
        "rt and i columns) or a matrix (with rt and i column)"
      )
    )
  }
  if (length(which(xic_mat$i > 0)) < minscan) {
    warning("Not enough point in xic to search for pics")
    return(NULL)
  }
  xic_mat <- xic_mat[order(rt)]
  ## remove trailing zero
  rt_range <- xic_mat[i > 0, range(rt)]
  xic_mat <- xic_mat[rt %between% rt_range]

  if (xic_mat[, .N] < minscan) {
    return(NULL)
  }

  xic_mat[, scan_nb := seq_len(.N)]

  peak_dt <- detect_centwave_peaks(
    vecint = xic_mat$i,
    majpeakonly = majpeakonly
  )

  if (isFALSE(peak_dt)) {
    return(NULL)
  }

  if (!is.null(minscan)) {
    peak_dt <- peak_dt[(scmax - scmin + 1) >= minscan, ]
  }

  peak_dt[, rtmin := xic_mat[scmin, rt]]
  peak_dt[, rtmax := xic_mat[scmax, rt]]
  peak_dt[, rt := xic_mat[scpos, rt]]
  peak_dt[, cpd_index := paste0("PK", seq_len(.N))]
  if (nrow(peak_dt) <= 0) {
    return(NULL)
  }
  output <- peak_dt[, {
    rt_range <- c(rtmin, rtmax)
    xic_i <- xic_mat[rt %between% rt_range]
    scan_nb <- xic_i[, .N]
    zigzag_score <- calculate_zigzag(xic_i[order(scan_nb), i])
    .(scan_nb = scan_nb, zigzag_score = as.numeric(zigzag_score))
  }, by = names(peak_dt)]

  return(output[])
}

#' Extract peaks from XICs
#'
#' @param xic_cpd a data.table with XIC informations:
#'                columns: `filter`, `EventIndex`, `rt`
#'                `i`
#' @inheritParams get_peaks_xic
#'
#' @export
#'
get_peaks_from_xic <- function(xic_cpd, minscan = 4) {
  if (!"CpdIndex" %in% names(xic_cpd)) {
    xic_cpd[, CpdIndex := 1]
  }
  if (!"filter" %in% names(xic_cpd) && "EventIndex" %in% names(xic_cpd)) {
    xic_cpd[, filter := EventIndex]
  }
  output <- xic_cpd[order(CpdIndex, rt), {
    peak_out <- get_peaks_xic(
      xic_mat = .SD[, .(rt, i)],
      minscan = minscan
    )
    peak_out
  }, by = .(CpdIndex, filter)]

  return(output)
}


#' Return best peak
#'
#' @param pk_scores A data.table containing peaks info and scores
#'
#' @return
#' A data.table containing the only best peaks for each `CpdIndex`
#' based on the maximum of `peak_score` column
#'
#' @export
#'
fun_best_peaks <- function(pk_scores) {
  if (!data.table::is.data.table(pk_scores)) {
    return(NULL)
  }
  pk_scores[, .SD[which.max(peak_score)], by = .(CpdIndex)]
}


#' Add peaks detection to annobject
#'
#' @param rt_limit limit rt deviation from reference when searching peaks
#' @inheritParams add_events
#' @inheritParams add_xics
#' @inheritParams get_peaks_xic
#'
#' @return
#' The annobject with the Peaks slots for each compound
#' filled.
#'
#' @export
#'
add_best_peaks <- function(
  annobject,
  minscan = 3,
  rt_limit = NULL,
  debug = FALSE
) {
  ## Select cpd with XICs
  cpd_to_get <- sapply(annobject$cpd, function(x) {
    nrow(x$XICs) > 0
  }) %>%
    which()
  for (i in cpd_to_get) {
    if (isTRUE(debug)) {
      message("cpd: ", i, " file: ", appendLF = FALSE)
    }
    best_peak_dt <- split(
      annobject$cpd[[i]]$XICs,
      annobject$cpd[[i]]$XICs$FileIndex
    ) %>%
      lapply(., function(x) {
        if (isTRUE(debug)) {
          message(x$FileIndex %>% unique(), "-", appendLF = FALSE)
        }
        ## Add CpdIndex in x
        CpdIndex_i <- i
        x[, CpdIndex := CpdIndex_i]
        xic_peaks <- get_peaks_from_xic(xic_cpd = x, minscan = minscan)
        if (!is.null(rt_limit) && nrow(xic_peaks) > 0) {
          rt_lim_range <- annobject$cpd[[i]]$cpd_info$rtmin +
            c(-rt_limit, +rt_limit)
          xic_peaks <- xic_peaks[rt %between% rt_lim_range]
        }

        pk_scores_dt <- get_peak_scores_mult(
          xic_peaks = xic_peaks,
          cpd_events = annobject$file$MSEvents[[unique(x$FileIndex)]][
            EventIndex %in%
              annobject$cpd[[i]]$MSEvents[
                FileIndex == unique(x$FileIndex), EventIndex
              ]
          ],
          CpdIndex = NULL,
          cpd_info = annobject$cpd[[i]]$cpd_info
        )
        if (nrow(pk_scores_dt) == 0) {
          return(NULL)
        } else {
          pk_best_dt <- fun_best_peaks(pk_scores_dt)
          pk_best_dt[, FileIndex := unique(x$FileIndex)]
          return(pk_best_dt[])
        }
      }) %>%
      rbindlist()
    if (isTRUE(debug)) {
      message("", appendLF = TRUE)
    }
    annobject$cpd[[i]]$Peaks <- best_peak_dt
  }
  return(annobject)
}
