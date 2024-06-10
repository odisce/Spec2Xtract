#' Save spectrum as .msp file
#'
#' @param spectra_dt Spectrum as data.table
#' @param spectra_info Spectrum info to add
#' @param cpd_info Compound info to add
#' @param file_out Output path `/path/to/file.msp`
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
      {. / 60} %>%
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
  info_to_add <- names(cpd_info)[!names(cpd_info) %in% c("CpdIndex", "compound", "mz", "rtsec", "inchikey", "elemcomposition", "rtmin", "mz_neutral", "mz_pos", "mz_neg")]
  if (length(info_to_add) > 0) {
    temp_add <- cpd_info[, ..info_to_add] %>%
      as.list()
    names(temp_add) <- stringr::str_to_upper(names(temp_add))
    info_list <- c(info_list[!names(info_list) %in% names(temp_add)], temp_add)
  }

  ## Add last 
  info_list$`NUM PEAKS` <- spectra_dt[, .N]

  ## Write msp
  file_out <- normalizePath(file_out)
  if (!dir.exists(dirname(file_out))) {dir.create(dirname(file_out))}
  if (file.exists(file_out)) {
    unlink(file_out)
  }
  ## Write header
  for (i in seq_len(length(info_list))) {
    str_to_append <- paste0(names(info_list[i]), ": ", info_list[[i]])
    output <- c(output, str_to_append)
    write(str_to_append, file = file_out, append = TRUE)
  }
  ## Write peak list
  spectra_dt <- spectra_dt[order(mz), ]
  for (i in seq_len(spectra_dt[, .N+2])) {
    if (i > spectra_dt[, .N]) {
      write("", file = file_out, append = TRUE)
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
      output <- c(output, str_to_append)
      write(str_to_append, file = file_out, append = TRUE)
    }
  }
  return(output)
}