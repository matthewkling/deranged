% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{landscape_template}
\alias{landscape_template}
\title{Generate example data structure for spatial data.}
\usage{
landscape_template(
  n_row = 10,
  n_col = 10,
  n_steps = 1,
  n_env = 1,
  names = c("s", "j", "a")
)
}
\arguments{
\item{n_row}{Number of rows in spatial grid.}

\item{n_col}{Number of columns in spatial grid.}

\item{n_steps}{Number of time steps for environmental variables.}

\item{n_env}{Number of environmental variables.}

\item{names}{Demographic stage names (character vector).}
}
\value{
A list of demographic and dispersal parameter objects.
}
\description{
Generate example data structure for spatial data.
}
