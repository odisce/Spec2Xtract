#' Extract a single XIC
#'
#' @inheritParams get_rawindex
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
