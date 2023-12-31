#' Check if raw file exist
#'
#' @param rawpath Path to the .raw file
#'
#' @return This function check .raw and .mzML files for readability.
#'         If the file is readable it will return the full path, if
#'         not it will stop with an error.
#'
#' @importFrom rawrr readTrailer
#' @importFrom tools file_ext
#' @importFrom MSnbase readMSData
#' @importFrom stringr str_to_lower
#' @export
#'
#' @examples
#' check_rawfile(rawrr::sampleFilePath())
check_rawfile <- function(rawpath) {
  if (!file.exists(rawpath)) {
    stop(paste0("Input file not found: ", rawpath))
  }
  switch(
    tools::file_ext(rawpath) %>% stringr::str_to_lower(),
    "raw" = {
      message("analyzing using .raw file implementation.")
      temp <- tryCatch({
        rawrr::readTrailer(rawpath)
      },
      error = function(e) {
        print(e)
        stop(
          paste0(
            "Cannot open the .raw file with rawrr. ",
            "The raw file must be from a Thermo MS spectrometer"
          )
        )
      })
      if (is.character(temp) && length(temp) == 1 && grepl("Unable to", temp)) {
        stop(
          paste0(
            "Cannot open the .raw file with rawrr. ",
            "The raw file must be from a Thermo MS spectrometer"
          )
        )
      }
      return(rawpath)
    },
    "mzml" = {
      message("analyzing using .mzML file implementation.")
      temp <- tryCatch(
        {
          MSnbase::readMSData(files = rawpath, mode = "onDisk")
        },
        error = function(e) {
          print(e)
          stop(
            paste0(
              "Cannot open the .mzML file with MSnbase. ",
              "Check compatibility or try a different format"
            )
          )
        }
      )
      return(rawpath)
    },
    "mzxml" = {
      message("analyzing using .mzXML file implementation.")
      temp <- tryCatch(
        MSnbase::readMSData(files = rawpath, mode = "onDisk"),
        error = function(e) {
          print(e)
          stop(
            paste0(
              "Cannot open the .mzXML file with MSnbase. ",
              "Check compatibility or try a different format"
            )
          )
        }
      )
      return(rawpath)
    }
  )
}
