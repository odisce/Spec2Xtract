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
  output <- index_table[, {
    msLevel <- gsub("^.*(ms)([0-9]{0,2}) .*$", "\\2", scanType)
    polarity <- gsub("^.*(\\+|\\-).*ms.*$", "\\1", scanType)
    spec_polarity <- ifelse(polarity == "+", 1, 0)
    if (msLevel == "") {
      spec_energy <- "NA"
      spec_coltype <- "NA"
      spec_prec <- "NA"
      msLevel <- 1
    } else {
      spec_prec <- gsub("^.*ms[0-9]{1,2} ([0-9]{1,5}\\.[0-9]{1,6})\\@.*$", "\\1", x = scanType)
      char_to_parse <- gsub("^.*\\@([A-z]{1,6}[0-9]{1,4}\\.[0-9]{0,4}).\\[.*$", "\\1", x = scanType)
      spec_energy <- gsub("([A-z]{1,6})([0-9]{1,6}\\.[0-9]{0,6})", "\\2", char_to_parse)
      spec_coltype <- gsub("([A-z]{1,6})([0-9]{1,6}\\.[0-9]{0,6})", "\\1", char_to_parse)
      msLevel <- as.numeric(msLevel)
    }
    .(
      msLevel = as.integer(msLevel),
      spec_energy = as.numeric(spec_energy),
      spec_coltype = as.character(spec_coltype),
      spec_polarity = as.integer(spec_polarity),
      spec_prec = as.numeric(spec_prec)
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
#' @export
#'
fun_check_cpd <- function(cpd) {
  x <- as.data.table(cpd)
  if (!any(c("compound", "rtsec", "elemcomposition") %in% names(x))) {
    stop("cpd_file format incorrect, need the following columns: compound, mz, rtsec, inchikey, elemcomposition")
  }
  if (!"inchikey" %in% names(x)) {
    warning("inchikey not found in cpd table, it may limit the export to other tools")
  }
  ## format file
  x[, cpd_iter := 1:.N]
  x[, compound := as.character(compound)]
  x[, rtsec := as.numeric(rtsec)]
  x[, rtmin := rtsec/60]
  x[, inchikey := as.character(inchikey)]
  x[, elemcomposition := as.character(elemcomposition) %>% stringr::str_trim()]
  return(x[])
}

#' Calculate ion mass from elemental composition
#'
#' @inheritParams fun_check_cpd
#' 
#' @import magrittr data.table
#' @importFrom Spec2Annot mz_from_string mz_calc_ion
#'
#' @return Return the \code{cpd} table with mz_neutral, mz_pos and mz_neg values calculated 
#'         from the elemental composition
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
  # files <- input_files
  # cpd <- cpd_dt
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
        return(F)
      }
    )
  })

  ## If any file is F, add info
  for (y in 1:nrow(temp_obj$file$info)) {
    if (!is.data.table(temp_obj$file$index[[y]]) && temp_obj$file$index[[y]] == F) {
      temp_obj$file$info[y, FileCheck := F]
    } else {
      temp_obj$file$info[y, FileCheck := T]
    }
  }

  ## Add compound list
  cpd <- fun_check_cpd(cpd)
  temp_obj$cpd <- lapply(1:nrow(cpd), function(x) {
    # x <- 1
    temp <- cpd[x,]
    cpd_dt_i <- tryCatch(
      {
        temp %>%
          cpd_add_ionsmass(.) %>% {
            .[, cpdCheck := T][]
          }
      }, error = function(e) {
        print(e)
        temp[, cpdCheck := F]
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