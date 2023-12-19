#' Extract a single XIC
#'
#' @inheritParams get_rawindex
#' @inheritParams fun_get_xics
#' @param filter special string to filter events
#'               (ex.: `"ms"`, or `"ms2"`)
#'
#' @return
#' A `data.table` containing the XIC with the following columns:
#'   `filter`: the scanType field
#'   `ppm`: the ppm parameter used
#'   `mass`: the center mass of the XIC
#'   `rt`: the retention time
#'   `i`: the intensities
#' @export
#' @examples
#' fun_get_xic(
#'   rawpath = rawrr::sampleFilePath(),
#'   filter = "ms2",
#'   mz = 669.8381,
#'   ppm = 10,
#'   type = "xic"
#' )
#'
fun_get_xic <- function(
  rawpath,
  mz,
  ppm,
  filter,
  type = c("xic", "tic", "bpc")[1]
) {
  output <- rawrr::readChromatogram(
    rawfile = rawpath,
    mass = mz,
    tol = ppm,
    filter = filter,
    type = type
  )
  if (!"times" %in% (names(output))) {
    output <- output[[1]]
  }
  out_dt <- data.table(
    "filter" = filter,
    "ppm" = as.numeric(output$ppm),
    "mass" = as.numeric(output$mass),
    "rt" = as.numeric(output$times),
    "i" = as.numeric(output$intensities)
  )
  return(out_dt)
}

#' Extract multiple XICs
#'
#' @inheritParams get_rawindex
#' @inheritParams get_events_types
#' @param mz centered mass of the XIC window
#' @param ppm mass range of the XIC window (in ppm)
#' @param mslevel MS(level) to extract as integer
#' @param isomz (optional) mass of the isolation window to extract
#' @param type string to extract either the `"xic"` or `"tic"`
#' @importFrom Spec2Annot mz_range
#'
#' @return
#' a data.table with the XICs
#' @export
#'
fun_get_xics <- function(
  rawpath,
  mz = NULL,
  ppm = 10,
  mslevel = 1L,
  isomz = NULL,
  index_table = NULL,
  type = c("xic", "tic")[1]
) {
  if (is.null(index_table)) {
    index_table <- get_rawindex(rawpath = rawpath)
  }
  if (mslevel == 1) {
    filter_vc <- index_table[MSOrder == "Ms", unique(scanType)]
  } else if (mslevel > 1) {
    if (!is.null(isomz)) {
      if (is.numeric(isomz)) {
        ## select only isomz events
        mz_range <- Spec2Annot::mz_range(isomz, ppm)
        filter_vc <- index_table[
          MSOrder == paste0("Ms", mslevel) &
            data.table::between(
              precursorMass,
              mz_range[1],
              mz_range[2]
            ),
          unique(scanType)
        ]
      } else if (is.character(isomz)) {
        if (isomz %in% index_table[, unique(scanType)]) {
          filter_vc <- isomz
        } else {
          stop(
            paste0(
              "isomz must be the m/Z value of the isolation window ",
              "or the exact name of the scanType"
            )
          )
        }
      }
    } else {
      filter_vc <- index_table[
        MSOrder == paste0("Ms", mslevel),
        unique(scanType)
      ]
    }
  }

  xics <- lapply(filter_vc, function(xx) {
    if (type == "tic") {
      output <- fun_get_xic(
        rawpath = rawpath,
        mz = NULL,
        ppm = NULL,
        filter = xx,
        type = type
      )
    } else {
      output <- fun_get_xic(
        rawpath = rawpath,
        mz = mz,
        ppm = ppm,
        filter = xx,
        type = type
      )
    }
    output[, type := type][, mslevel := mslevel]
    return(output[])
  }) %>%
    rbindlist()
  return(xics)
}

#' Add XICs for all MS(n) events and compounds
#'
#' @inheritParams add_events
#' @param ms1 Logical
#' @param debug Logical to show log messages (`TRUE`) or not (`FALSE`)
#'
#' @return
#' Populate the `XICs` fields for each compound
#' in the annobject using information from the cpd
#' table.
#'
#' @export
#'
add_xics <- function(
  annobject,
  ms1 = TRUE,
  debug = FALSE
) {
  for (i in seq_len(length(annobject$cpd))) {
    if (isTRUE(debug)) {
      message("[XICs] CPD", i, " ", appendLF = FALSE)
    }

    all_xics <- lapply(
      seq_len(nrow(annobject$file$info)),
      function(file_i) {
        if (isTRUE(debug)) {
          message("F", file_i, appendLF = FALSE)
        }
        if (
          annobject$file$info[file_i, FileCheck] == FALSE |
            (
              is.null(
                annobject$cpd[[i]]$MSEvents[FileIndex == file_i]
              ) ||
                nrow(annobject$cpd[[i]]$MSEvents[FileIndex == file_i]) == 0
            )
        ) {
          next
        }
        cpd_events_i <- merge(
          annobject$cpd[[i]]$MSEvents[FileIndex == file_i],
          annobject$file$MSEvents[[file_i]],
          by = "MSEvent_index"
        )
        if (isTRUE(ms1) & cpd_events_i[, any(msLevel == 1)]) {
          cpd_events_i <- cpd_events_i[msLevel == 1, ]
        }
        ## For each events
        events_max <- nrow(cpd_events_i)
        file_xics <- lapply(
          seq_len(events_max),
          function(events_i) {
            if (isTRUE(debug)) {
              message("E", events_i, appendLF = FALSE)
              if (events_i < events_max) {
                message("-", appendLF = FALSE)
              }
            }
            cpd_events_ii <- cpd_events_i[events_i, ]
            temp_xic <- fun_get_xics(
              rawpath = annobject$file$info[file_i, file_path],
              mz = cpd_events_ii[,
                ifelse(
                  spec_polarity == 0,
                  annobject$cpd[[i]]$cpd_info$mz_neg,
                  annobject$cpd[[i]]$cpd_info$mz_pos
                )
              ],
              ppm = 10,
              mslevel = cpd_events_ii[, msLevel],
              isomz = cpd_events_ii[, scanType],
              index_table = annobject$file$index[[file_i]],
              type = "xic"
            )
            temp_xic[, MSEvent_index := cpd_events_ii[, MSEvent_index]]
            return(temp_xic[, .(ppm, rt, i, type, MSEvent_index)])
          }
        ) %>%
          rbindlist()
        file_xics[, FileIndex := file_i]
        return(file_xics)
      }
    ) %>%
      rbindlist()
    annobject$cpd[[i]]$XICs <- all_xics
    if (isTRUE(debug)) {
      message("", appendLF = TRUE)
    }
  }
  return(annobject)
}
