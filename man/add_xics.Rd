% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_xics.R
\name{add_xics}
\alias{add_xics}
\title{Add XICs for all MS(n) events and compounds}
\usage{
add_xics(annobject, ms1 = TRUE, debug = FALSE)
}
\arguments{
\item{annobject}{Object build with \code{Spec2Extract::init_object()}}

\item{ms1}{Logical}

\item{debug}{Logical to show log messages (\code{TRUE}) or not (\code{FALSE})}
}
\value{
Populate the \code{XICs} fields for each compound
in the annobject using information from the cpd
table.
}
\description{
Add XICs for all MS(n) events and compounds
}
