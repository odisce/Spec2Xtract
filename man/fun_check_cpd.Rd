% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_object.R
\name{fun_check_cpd}
\alias{fun_check_cpd}
\title{Check compound file}
\usage{
fun_check_cpd(cpd)
}
\arguments{
\item{cpd}{A table containing compound information, must have the following
columns: compound, rtsec, elemcomposition and optionally inchikey}
}
\value{
Return a formated data.table with compound informations
}
\description{
Check compound file
}
\examples{
fun_check_cpd(Spec2Xtract:::example_cpdlist_realdt)

}
