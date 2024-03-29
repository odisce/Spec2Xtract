% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_best_peaks.R
\name{get_peak_scores}
\alias{get_peak_scores}
\title{Calculate peak scores}
\usage{
get_peak_scores(peaks_dt, ref_rtmin = NULL)
}
\arguments{
\item{peaks_dt}{Peak list as returned by
\code{Spec2Xtract::get_peaks_from_xic()}}

\item{ref_rtmin}{reference rt in minute}
}
\value{
A data.table identical to \code{peaks_dt} with the following
supplementary columns:
\itemize{
\item \code{peak_in_range}: Number of peaks in range for this MS Event
\item \code{diff_from_ref}: RT shift from reference
\item \code{peak_score}: Final score
}
}
\description{
Calculate peak scores
}
