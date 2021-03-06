% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_mods.R
\name{pkmod1cpt}
\alias{pkmod1cpt}
\title{One compartment IV infusion with first-order elimination.}
\usage{
pkmod1cpt(tm, kR, pars, init = 0, inittm = 0)
}
\arguments{
\item{tm}{Vector of times to evaluate the PK function at}

\item{kR}{Infusion rate (e.g. ml/min).}

\item{pars}{Named vector of parameters with names ('ke','v') or ('cl').}

\item{init}{Initial concentration. Defaults to 0.}

\item{inittm}{Time of initiation of infusion. Defaults to 0.}
}
\value{
Returns a vector of numeric values describing predicted concentrations
for a one-compartment model with IV infusion.
}
\description{
One compartment IV infusion with first-order elimination.
}
\examples{
pkmod1cpt(1,1,c(ke = 0.5, v = 1))
}
