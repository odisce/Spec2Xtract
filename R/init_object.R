#' Get .raw index as data.table
#'
#' @inheritParams check_rawfile
#'
#' @importFrom rawrr readIndex
#' @import magrittr data.table
#'
#' @return
#' Return a data.table with the spectra indexes from \code{rawrr::readIndex}.
#'
#' @export
#'
#' @examples
#' get_rawindex(rawrr::sampleFilePath())
get_rawindex <- function(rawpath) {
  check_rawfile(rawpath)
  rawrr::readIndex(rawfile = rawpath) %>%
    as.data.table() %>%
    return()
}

#' Extract ms level and collision type, energy from index
#'
#' @param index_table Index table returned by \code{get_rawindex}
#'
#' @import data.table magrittr
#'
#' @return
#' Add the following columns to \code{index_table}:
#'   - `"msLevel"`: Fragmentation level
#'   - `"spec_energy"`: Collision energy
#'   - `"spec_coltype"`: Collision type
#'   - `"spec_polarity"`: Polarity
#'   - `"spec_prec"`: Precursor m/Z
#'
#' @export
#'
parse_index_dt <- function(index_table) {

  pattern_search <- list(
    "spec_prec" = list(
      "pattern" = "^.*ms[0-9]{1,2} ([0-9]{1,5}\\.[0-9]{1,6})\\@.*$",
      "subset" = "\\1"
    ),
    "char_to_parse" = list(
      "pattern" = "^.*\\@([a-zA-Z]{1,6}[0-9]{1,4}\\.[0-9]{0,4}).\\[.*$",
      "subset" = "\\1"
    ),
    "spec_energy" = list(
      "pattern" = "^.*\\@([a-zA-Z]{1,6})([0-9]{1,6}\\.[0-9]{0,6}).\\[.*$",
      "subset" = "\\2"
    ),
    "spec_coltype" = list(
      "pattern" = "^.*\\@([a-zA-Z]{1,6})([0-9]{1,6}\\.[0-9]{0,6}).\\[.*$",
      "subset" = "\\1"
    )
  )

  output <- index_table[, {
    msLevel <- gsub("^.*(ms)([0-9]{0,2}) .*$", "\\2", scanType)
    if (msLevel == "") {
      pattern_res <- list(
        "spec_energy" = FALSE,
        "spec_coltype" = FALSE,
        "spec_prec" = FALSE
      )
      pattern_res$mslevel <- 1
    } else {
      pattern_res <- lapply(pattern_search, function(x) {
        if (grepl(pattern = x$pattern, x = scanType)) {
          gsub(pattern = x$pattern, replacement = x$subset, x = scanType)
        } else {
          return(FALSE)
        }
      })
      pattern_res$mslevel <- msLevel
      polarity <- gsub("^.*(\\+|\\-).*ms.*$", "\\1", scanType)
      pattern_res$spec_polarity <- ifelse(polarity == "+", 1, 0)
    }
    .(
      msLevel = as.integer(pattern_res$mslevel),
      spec_energy = pattern_res$spec_energy %>% {
        ifelse(isFALSE(.), as.numeric(NA), as.numeric(.))
      },
      spec_coltype = pattern_res$spec_coltype %>% {
        ifelse(isFALSE(.), as.character(NA), as.character(.))
      },
      spec_polarity = pattern_res$spec_polarity %>% {
        ifelse(isFALSE(.), as.integer(NA), as.integer(.))
      },
      spec_prec = pattern_res$spec_prec %>% {
        ifelse(isFALSE(.), as.numeric(NA), as.numeric(.))
      }
    )
  }, by = .(scan)]
  merge(index_table, output, by = "scan") %>%
    return()
}

#' Check compound file
#'
#' @param cpd A table containing compound information, must have the following
#'            columns: compound, rtsec, elemcomposition and optionally inchikey
#'
#' @import data.table magrittr
#' @return Return a formated data.table with compound informations
#' @examples
#' fun_check_cpd(Spec2Xtract:::example_cpdlist)
#'
#' @export
#'
fun_check_cpd <- function(cpd) {
  x <- as.data.table(cpd)
  if (!any(c("compound", "rtsec", "elemcomposition") %in% names(x))) {
    stop(
      paste0(
        "cpd_file format incorrect, need the following columns: ",
        "compound, mz, rtsec, inchikey, elemcomposition"
      )
    )
  }
  if (!"inchikey" %in% names(x)) {
    warning(
      paste0(
        "inchikey not found in cpd table, ",
        "it may limit the export to other tools"
      )
    )
    inchikey <- as.character(NA)
  }
  ## format file
  x[, cpd_iter := seq_len(.N)]
  x[, compound := as.character(compound)]
  x[, rtsec := as.numeric(rtsec)]
  x[, rtmin := rtsec / 60]
  x[, inchikey := as.character(inchikey)]
  x[, elemcomposition := as.character(elemcomposition) %>% stringr::str_trim()]
  return(x[])
}

#' Calculate ion mass from elemental composition
#'
#' @inheritParams fun_check_cpd
#'
#' @import magrittr data.table Spec2Annot
#' @importFrom Spec2Annot mz_from_string mz_calc_ion
#'
#' @return Return the \code{cpd} table with:
#'         mz_neutral, mz_pos and mz_neg values
#'         calculated from the elemental composition
#' @export
#'
cpd_add_ionsmass <- function(cpd) {
  temp <- cpd[, {
    mz_neutral <- Spec2Annot::mz_from_string(string = elemcomposition)
    mz_neg <- Spec2Annot::mz_calc_ion(mz_neutral, form = "-H")
    mz_pos <- Spec2Annot::mz_calc_ion(mz_neutral, form = "+H")
    .(
      mz_neutral = mz_neutral,
      mz_pos = mz_pos,
      mz_neg = mz_neg
    )
  }, by = cpd_iter]
  merge(cpd, temp, by = "cpd_iter") %>%
    return()
}

#' Initialize object to store results
#'
#' @param files Vector containing the paths to the .raw or .mzML files
#' @inheritParams fun_check_cpd
#'
#' @import data.table magrittr
#'
#' @return
#' Return a list with the following levels:
#'   - `"cpd_info"`: data.table with information on the compounds from cpd.
#'   - `"XICs"`: level to store XICs informations
#'   - `"Peaks"`: level to store peaks informations
#'   - `"MSspectra"`: level to store mass spectrum
#'
#' @export
#'
init_object <- function(files, cpd) {
  ## Add files info
  temp_obj <- list(
    "file" = list(
      "info" = data.table(
        file_path = files,
        file_name = files %>%
          basename(.) %>%
          tools::file_path_sans_ext(.)
      )
    )
  )
  ## Add files index
  temp_obj$file$index <- lapply(files, function(x) {
    raw_index <- tryCatch(
      {
        x %>%
          get_rawindex() %>%
          parse_index_dt()
      }, error = function(e) {
        print(e)
        return(FALSE)
      }
    )
    return(raw_index)
  })

  ## If any file is F, add info
  for (y in seq_len(nrow(temp_obj$file$info))) {
    if (!is.data.table(temp_obj$file$index[[y]]) &&
          temp_obj$file$index[[y]] == FALSE) {
      temp_obj$file$info[y, FileCheck := FALSE]
    } else {
      temp_obj$file$info[y, FileCheck := TRUE]
    }
  }

  ## Add compound list
  cpd <- fun_check_cpd(cpd)
  temp_obj$cpd <- lapply(seq_len(nrow(cpd)), function(x) {
    temp <- cpd[x, ]
    cpd_dt_i <- tryCatch({
      temp %>%
        cpd_add_ionsmass(.) %>%
        {
          .[, cpdCheck := TRUE][]
        }
    }, error = function(e) {
      print(e)
      temp[, cpdCheck := FALSE]
      temp[]
    })

    list(
      "cpd_info" = cpd_dt_i,
      "XICs" = NULL,
      "Peaks" = NULL,
      "MSspectra" = NULL
    )
  })

  return(temp_obj)
}
