#' Write file to disk
#'
#' @param string string as vector or list to write to the file
#' @param file_path Path to the output file
#' @param overwrite Boolean to append `FALSE` or overwrite
#' the file `TRUE`
#'
#' @return
#' A Boolean if the file was written `TRUE` or not `FALSE`
#'
#' @export
write_file <- function(string, file_path, overwrite = TRUE) {
  file_out <- normalizePath(file_path, mustWork = FALSE)
  if (!dir.exists(dirname(file_out))) {
    dir.create(dirname(file_out))
  }
  if (isTRUE(overwrite) && file.exists(file_out)) {
    unlink(file_out)
  }
  if (is.list(string)) {
    string <- unlist(string)
  }
  write(string, file = file_out, append = TRUE)
  if (file.exists(file_out)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Save spectrum as .msp file
#'
#' @param spectra_dt Spectrum as data.table with columns
#' (mz, irel, formula)
#' @param spectra_info Spectrum info to add
#' @param cpd_info Compound info to add
#' @param file_out Output path (example: `/path/to/file.msp`)
#' if `NULL` the function will return a list with the file
#' content.
#'
#' @importFrom stringr str_to_upper
#'
#' @return
#' Save the .msp file to the file_out directory and
#' return a list with the .msp content
#'
#' @export
#'
save_as_msp <- function(
  spectra_dt,
  spectra_info,
  cpd_info,
  file_out
) {
  ## Checks
  if (!is.data.table(spectra_dt) && !all(c("mz", "irel", "formula") %in% names(spectra_dt))) {
    stop("save_as_msp: spectra_dt input is not conform")
  }
  if (
    !is.data.table(spectra_info) &&
      !all(c(
        "precursorMass",
        "spec_polarity",
        "spec_rt",
        "spec_polarity",
        "spec_energy",
        "spec_coltype",
        "msLevel",
        "spec_isowin",
        "isopurity"
      ) %in% names(spectra_info))
  ) {
    stop("save_as_msp: spectra_info input is not conform")
  }

  if (
    !is.data.table(cpd_info) &&
      !all(c(
        "compound",
        "elemcomposition",
        "inchikey",
        "mz_neutral"
      ) %in% names(spectra_dt))
  ) {
    stop("save_as_msp: cpd_info input is not conform")
  }
  output <- c()
  info_list <- list(
    "NAME" = cpd_info$compound %>% as.character(),
    "PRECURSORMZ" = spectra_info$precursorMass %>% sprintf("%0.6f", .),
    "PRECURSORTYPE" = if (!"PRECURSORTYPE" %in% names(cpd_info)) {
      if (spectra_info$spec_polarity == 1) {
        "[M+H]+"
      } else {
        "[M-H]-"
      }
    },
    "FORMULA" = cpd_info$elemcomposition %>% as.character(),
    "INCHIKEY" = cpd_info$inchikey %>% as.character(),
    "RETENTIONTIME" = spectra_info$spec_rt %>%
      as.numeric() %>%
      {
        . / 60
      } %>%
      sprintf("%.2f", .),
    "IONMODE" = ifelse(spectra_info$spec_polarity == 1, "positive", "negative"),
    "COLLISIONENERGY" = spectra_info$spec_energy,
    "FRAGMENTATIONMODE" = spectra_info$spec_coltype %>%
      stringr::str_to_upper(),
    "EXACTMASS" = cpd_info$mz_neutral,
    "MSLEVEL" = spectra_info$msLevel,
    "ISOLATIONWIN" = spectra_info$spec_isowin,
    "COMMENT" = paste0(
      "Extracted by Spec2Xtract",
      " ; ",
      "IsoPurity: ",
      sprintf("%0.2f", as.numeric(spectra_info$isopurity))
    )
  )

  ## Replace by cpd_info column names if present
  info_to_add <- names(cpd_info)[
    !names(cpd_info) %in% c(
      "CpdIndex",
      "compound",
      "mz",
      "rtsec",
      "inchikey",
      "elemcomposition",
      "rtmin",
      "mz_neutral",
      "mz_pos",
      "mz_neg"
    )
  ]
  if (length(info_to_add) > 0) {
    temp_add <- cpd_info[, ..info_to_add] %>%
      as.list()
    names(temp_add) <- stringr::str_to_upper(names(temp_add))
    info_list <- c(info_list[!names(info_list) %in% names(temp_add)], temp_add)
  }

  ## Add last
  info_list$`NUM PEAKS` <- spectra_dt[, .N]

  ## Generate content
  ## Write header
  for (i in seq_len(length(info_list))) {
    str_to_append <- paste0(names(info_list[i]), ": ", info_list[[i]])
    output <- c(output, str_to_append)
  }
  ## Write peak list
  spectra_dt <- spectra_dt[order(mz), ]
  for (i in seq_len(spectra_dt[, .N + 2])) {
    if (i > spectra_dt[, .N]) {
      str_to_append <- ""
    } else {
      str_to_append <- paste0(
        spectra_dt[i, mz] %>% sprintf("%.6f", .), " ",
        spectra_dt[i, irel] %>% sprintf("%.4f", .)
      )
      if ("formula" %in% names(spectra_dt) && spectra_dt[i, !is.na(formula)]) {
        str_to_append <- paste0(
          str_to_append, " ",
          spectra_dt[i, formula] %>% as.character()
        )
      }
    }
    output <- c(output, str_to_append)
  }
  if (!is.null(file_out)) {
    ## Write msp
    write_file(
      string = output,
      file_path = file_out,
      overwrite = TRUE
    )
  }
  return(output)
}
