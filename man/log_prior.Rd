% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{log_prior}
\alias{log_prior}
\title{Calculate logged-prior probability for a set of parameters}
\usage{
log_prior(lpars, pkmod)
}
\arguments{
\item{lpars}{Logged parameter values. Can be a subset of the full set of PK or PK-PD parameter values.}

\item{pkmod}{`pkmod` object. Mean values are a subset of log(pars_pk), log(pars_pd),
log(sigma_add), log(sigma_mult). PK-PD parameter values not specified in `lpars` will be inferred from `pkmod`.}
}
\value{
Numeric value of length 1
}
\description{
Calculate logged-prior probability for a set of parameters, assuming that parameter values
are log-normally distributed. Mean values are set as the logged parameter values in
the `pkmod` object. Variances are given by the diagonal elements of `prior_vcov`.
}
\examples{
# evaluate log-prior for pk parameters + residual
lpars = log(c(cl=11,q2=3,q3=25,v=15,v2=30,v3=50,ke0=1.15,sigma_add=0.15))
prior_vcov <- matrix(diag(c(0.265,0.346,0.209,0.610,0.565,0.597,0.702,0.463)), 8,8,
dimnames = list(NULL,names(lpars)))
my_mod <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2),
sigma_add = 0.2, log_response = TRUE, Omega = prior_vcov)
log_prior(lpars, my_mod)

# evaluate log-prior for subset of parameters (remove volume parameters)
lpars_sub = log(c(cl=11,q2=3,q3=25,ke0=1.15,sigma_add=0.15))
prior_vcov_sub <- matrix(diag(c(0.265,0.346,0.209,0.702,0.463)), 5,5,
dimnames = list(NULL,names(lpars_sub)))
my_mod <- update(my_mod, Omega = prior_vcov_sub)
log_prior(lpars_sub, my_mod)
}
