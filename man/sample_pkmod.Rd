% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poppk_mods.R
\name{sample_pkmod}
\alias{sample_pkmod}
\title{Sample parameters from a `pkmod` object}
\usage{
sample_pkmod(pkmod, log_normal = TRUE, ...)
}
\arguments{
\item{pkmod}{`pkmod` object with associated Omega matrix describing random effect variances}

\item{log_normal}{Logical. Assumes random effects are log-normally distributed
and multiplicative if TRUE, additive and normally distributed if FALSE.}

\item{...}{Arguments passed to update.pkmod}
}
\value{
`pkmod` object with updated parameters.
}
\description{
Sample parameters from a `pkmod` object
}
\examples{
sample_pkmod(pkmod_schnider(AGE = 40,HGT=170,TBW=50,MALE=TRUE))
sample_pkmod(pkmod_eleveld_ppf(AGE = 40,TBW = 56,HGT=150,MALE = TRUE, PD = FALSE))
}
