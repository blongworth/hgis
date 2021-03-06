% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hgis_data_reduction_functions.R
\name{reduce_hgis}
\alias{reduce_hgis}
\title{Reduce HGIS data from an AMS results file}
\usage{
reduce_hgis(
  file,
  date = NULL,
  standards = NULL,
  blanks = NULL,
  outliers = NULL,
  remove_outliers = TRUE,
  get_consensus = TRUE
)
}
\arguments{
\item{file}{A SNICS results file}

\item{date}{Date sample run if analysing a specific day}

\item{standards}{A vector of standard positions.}

\item{blanks}{A vector of blank positions.}
}
\value{
A dataframe of results.
}
\description{
Steps:
\enumerate{
\item Load data
\item Process raw data.
\item Insert raw data into DB
\item Reduce data
\item Normalize
\item Blank Correct
\item Insert results into DB
}
}
