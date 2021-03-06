% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pkpd_assist.R
\name{eleveld_vcov}
\alias{eleveld_vcov}
\title{Generate variance-covariance matrix for Eleveld PK-PD model}
\usage{
eleveld_vcov(
  dat,
  N = 1000,
  rates = TRUE,
  varnames = c("K10", "K12", "K21", "K13", "K31", "V1", "V2", "V3", "KE0", "CE50",
    "SIGMA")
)
}
\arguments{
\item{dat}{Data frame of observed patient covariates}

\item{N}{Number of Monte Carlo samples}

\item{rates}{Logical. Should rate constants be calculated}

\item{varnames}{Column names of variables used to calculate variance-covariance matrix}
}
\value{
Returns a variance-covariance matrix as a numeric matrix.
}
\description{
Generate the variance-covariance matrix for Eleveld PK-PD model for an observation
via Monte Carlo sampling.
}
