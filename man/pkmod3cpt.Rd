% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_mods.R
\name{pkmod3cpt}
\alias{pkmod3cpt}
\title{Three compartment IV infusion with first-order elimination.}
\usage{
pkmod3cpt(tm, kR, pars, init = c(0, 0, 0), inittm = 0, k20 = 0, k30 = 0)
}
\arguments{
\item{tm}{Vector of times to evaluate the PK function at}

\item{kR}{Infusion rate (e.g. ml/min).}

\item{pars}{Named vector of parameters with names ('K10','K12','K21','V1','V2') or ('CL','Q','V1','V2').}

\item{init}{Initial concentration. Defaults to 0 in all compartments.}

\item{inittm}{Time of initiation of infusion. Defaults to 0.}

\item{k20}{Elimination rate constant for second compartment. Defaults to 0.}

\item{k30}{Elimination rate constant for second compartment. Defaults to 0.}
}
\value{
Returns a numeric matrix with predicted concentrations
for a three-compartment model with IV infusion. Compartments are represented in
rows while times are in columns.
}
\description{
Three compartment IV infusion with first-order elimination.
}
\examples{
pkmod3cpt(1,1,c(CL = 15, Q2 = 10, Q3 = 5, V1 = 10, V2 = 20, V3 = 50))
}
