# --------------------------------------------------------------------------------------------------------------------------------
# - Library of PK functions for 1, 2, and 3 compartment models with effect-site compartment for 3 cmpt model ---------------------
# --------------------------------------------------------------------------------------------------------------------------------

#' One compartment IV infusion with first-order elimination.
#'
#'
#'
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names ('ke','v') or ('cl').
#' @param init Initial concentration. Defaults to 0.
#' @param inittm Time of initiation of infusion. Defaults to 0.
#' @return Returns a vector of numeric values describing predicted concentrations
#' for a one-compartment model with IV infusion.
#'
#' @examples
#' pkmod1cpt(1,1,c(ke = 0.5, v = 1))
#' @export
pkmod1cpt <- function(tm, kR, pars, init = 0, inittm = 0){

  names(pars) <- tolower(names(pars))
  tm <- tm - inittm

  if(any(!(c("ke","v") %in% names(pars))) & any(!(c("cl","v") %in% names(pars))))
    stop("pars must have names ('ke','v') or ('cl','v')")

  v <- pars["v"]
  if(!("ke" %in% names(pars)))
    ke <- pars["cl"] / v
  else
    ke <- pars["ke"]

  return((kR/ke*(1-exp(-tm*ke)) + init*v * exp(-tm*ke)) / v)
}
class(pkmod1cpt) <- "pkmod"




#' Two compartment IV infusion with first-order elimination.
#'
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names ('K10','K12','K21','V1','V2') or ('CL','Q','V1','V2').
#' @param init Initial concentration. Defaults to 0 in both compartments.
#' @param inittm Time of initiation of infusion. Defaults to 0.
#' @param k20 Elimination rate constant for second compartment. Defaults to 0.
#' @return Returns a numeric matrix with predicted concentrations
#' for a two-compartment model with IV infusion. Compartments are represented in
#' rows while times are in columns.
#' @examples
#' pkmod2cpt(1,1,c(CL = 15, V1 = 10, Q = 10, V2 = 20))
#' @export
pkmod2cpt <- function(tm, kR, pars, init = c(0,0), inittm = 0, k20 = 0){

  names(pars) <- tolower(names(pars))
  tm <- tm - inittm

  if(any(!(c("k10","k12","k21","v1","v2") %in% names(pars))) & any(!(c("cl","q","v1","v2") %in% names(pars))))
    stop("pars must have names ('k10','k12','k21','v1','v2') or ('cl','q','v1','v2')")

  v1 <- pars["v1"]
  v2 <- pars["v2"]
  if(all(c("cl","q") %in% names(pars))){
    k10 <- pars["cl"] / v1
    k12 <- pars["q"] / v1
    k21 <- pars["q"] / v2
  } else{
    k10 <- pars["k10"]
    k12 <- pars["k12"]
    k21 <- pars["k21"]
  }

  E1 <- k10+k12
  E2 <- k21+k20

  lambda1 = 0.5*((E1+E2)+sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))
  lambda2 = 0.5*((E1+E2)-sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  Doserate <- kR

  A1term1 = (((A1last*E2+Doserate+A2last*k21)-A1last*lambda1)*exp(-tm*lambda1)-((A1last*E2+Doserate+A2last*k21)-A1last*lambda2)*exp(-tm*lambda2))/(lambda2-lambda1)
  A1term2 = Doserate*E2*(1/(lambda1*lambda2)+exp(-tm*lambda1)/(lambda1*(lambda1-lambda2))-exp(-tm*lambda2)/(lambda2*(lambda1-lambda2)))
  A1term <- A1term1+A1term2 # Amount in central compartment

  A2term1 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-tm*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-tm*lambda2))/(lambda2-lambda1)
  A2term2 = Doserate*k12*(1/(lambda1*lambda2)+exp(-tm*lambda1)/(lambda1*(lambda1-lambda2))-exp(-tm*lambda2)/(lambda2*(lambda1-lambda2)))
  A2term <- A2term1+A2term2   #Amount in the peripheral compartment

  return(rbind(A1term/v1, A2term/v2))
}
class(pkmod2cpt) <- "pkmod"




#' Three compartment IV infusion with first-order elimination.
#'
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names ('K10','K12','K21','V1','V2') or ('CL','Q','V1','V2').
#' @param init Initial concentration. Defaults to 0 in all compartments.
#' @param inittm Time of initiation of infusion. Defaults to 0.
#' @param k20 Elimination rate constant for second compartment. Defaults to 0.
#' @param k30 Elimination rate constant for second compartment. Defaults to 0.
#' @return Returns a numeric matrix with predicted concentrations
#' for a three-compartment model with IV infusion. Compartments are represented in
#' rows while times are in columns.
#' @examples
#' pkmod3cpt(1,1,c(CL = 15, Q2 = 10, Q3 = 5, V1 = 10, V2 = 20, V3 = 50))
#' @export
pkmod3cpt <- function(tm, kR, pars, init = c(0,0,0), inittm = 0, k20 = 0, k30 = 0){

  names(pars) <- tolower(names(pars))
  tm <- tm - inittm

  if(any(!(c("k10","k12","k21","k13","k31","v1","v2","v3") %in% names(pars))) & any(!(c("cl","q2","q3","v1","v2","v3") %in% names(pars))))
    stop("pars must have names ('k10','k12','k21','k13','k31','v1','v2','v3') or ('cl','q2','q3','v1','v2','v3')")

  v1 <- pars["v1"]
  v2 <- pars["v2"]
  v3 <- pars["v3"]

  if(all(c("cl","q2","q3") %in% names(pars))){
    k10 <- pars["cl"] / v1
    k12 <- pars["q2"] / v1
    k21 <- pars["q2"] / v2
    k13 <- pars["q3"] / v1
    k31 <- pars["q3"] / v3
  } else{
    k10 <- pars["k10"]
    k12 <- pars["k12"]
    k21 <- pars["k21"]
    k13 <- pars["k13"]
    k31 <- pars["k31"]
  }

  E1 <- k10+k12+k13
  E2 <- k21+k20
  E3 <- k31+k30

  #calculate hybrid rate constants
  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(pmax(0,-1*Q))
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  A3last <- init[3]*v3
  Doserate <- kR

  B = A2last*k21+A3last*k31
  C = E3*A2last*k21+E2*A3last*k31
  I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
  J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

  A1term1 = A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
  A1term2 = exp(-tm*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
  A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
  A1term  = A1term1+A1term2+A1term3    #Amount in the central compartment

  A2term1 = A2last*(exp(-tm*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
  A2term2 = exp(-tm*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
  A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
  A2term  = A2term1+A2term2+A2term3    #Amount in the first-peripheral compartment

  A3term1 = A3last*(exp(-tm*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
  A3term2 = exp(-tm*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
  A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
  A3term = A3term1+A3term2+A3term3

  return(rbind(A1term/v1, A2term/v2, A3term/v3))
}
class(pkmod3cpt) <- "pkmod"




#' Solution to three-compartment IV model with effect-site
#'
#' 3 compartment IV infusion with first-order absorption between compartments and with an additional effect-site compartment.
#' The analytical solutions implemented in this function are provided in "ADVAN-style analytical solutions for common pharmacokinetic models" by
#' Abuhelwa et al. 2015.
#'
#' This function takes in arguments for each of the absorption and elimination rate constants of a three-compartment model
#' as well as initial concentrations, c0. ke0 gives the rate of elimination from the effect-site compartment into the
#' central compartment (i.e. k41). The rate of absorption into the effect-site compartment is set at 1/10,000 the value of ke0.
#' The function returns a set of functions that calculate the concentration in each of the four compartments as a function of
#' time.
#' @param tm Vector of times to evaluate the PK function at
#' @param kR Infusion rate (e.g. ml/min).
#' @param pars Named vector of parameters with names (k10,k12,k21,k13,k31,v1,v2,v3,ke0)
#' @param init Initial concentration
#' @param inittm Time of initiation of infusion
#' @param returncpt Optionally specify a single compartment to return concentrations for.
#' Defaults to returning all compartment concentrations.
#' @return Returns a numeric matrix with predicted concentrations
#' for a three-compartment model with a metabolite/effect-site compartment and IV infusion.
#' Compartments are represented in rows while times are in columns.
#' @examples
#' pars_3cpt <- c(k10=1.5,k12=0.15,k21=0.09,k13=0.8,k31=0.8,v1=10,v2=15,v3=100,ke0=1)
#' pkmod3cptm(1,1,pars_3cpt)
#' @export
pkmod3cptm <- function(tm, kR, pars, init = c(0,0,0,0), inittm = 0,
                       returncpt = c("all","cpt1","cpt2","cpt3","cpt4")) {

  names(pars) <- tolower(names(pars))

  if(any(!(c("k10","k12","k21","k13","k31","ke0","v1","v2","v3") %in% names(pars))) &
     any(!(c("cl","q2","q3","ke0","v1","v2","v3") %in% names(pars))))
    stop("pars must have names ('k10','k12','k21','k13','k31','v1','v2','v3','ke0') or ('cl','q2','q3','v1','v2','v3','ke0')")

  v1 <- pars["v1"]
  v2 <- pars["v2"]
  v3 <- pars["v3"]
  kme <- ke0 <- pars["ke0"]

  if(all(c("cl","q2","q3") %in% names(pars))){
    k10 <- pars["cl"] / v1
    k12 <- pars["q2"] / v1
    k21 <- pars["q2"] / v2
    k13 <- pars["q3"] / v1
    k31 <- pars["q3"] / v3
  } else{
    k10 <- pars["k10"]
    k12 <- pars["k12"]
    k21 <- pars["k21"]
    k13 <- pars["k13"]
    k31 <- pars["k31"]
  }

  if(!("k20" %in% names(pars))){
    k20 <- 0
  } else{
    k20 <- pars["k20"]
  }
  if(!("k30" %in% names(pars))){
    k30 <- 0
  } else{
    k30 <- pars["k30"]
  }
  if(!("km" %in% names(pars))){
    km  <- kme / 1e5
  } else{
    km <- pars["km"]
  }
  if(!("v4" %in% names(pars))){
    v4  <- v1 / 1e5
  } else{
    v4 <- pars["v4"]
  }

  returncpt <- match.arg(returncpt)
  tm <- tm - inittm

  E1 <- k10+k12+k13+km
  E2 <- k21+k20
  E3 <- k31+k30

  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(pmax(0,-1*Q))
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

  A1last <- init[1]*v1
  A2last <- init[2]*v2
  A3last <- init[3]*v3
  Amlast <- init[4]*v4

  B = A2last*k21+A3last*k31
  C = E3*A2last*k21+E2*A3last*k31
  I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
  J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

  if(returncpt %in% c("all", "cpt1")){
    A1term1 = A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-tm*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A1term3 = kR*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A1term = A1term1+A1term2+A1term3
  } else A1term = NULL


  if(returncpt %in% c("all", "cpt2")){
    A2term1 = A2last*(exp(-tm*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-tm*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = kR*k12*(E3/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A2term = A2term1+A2term2+A2term3
  } else A2term = NULL


  if(returncpt %in% c("all", "cpt3")){
    A3term1 = A3last*(exp(-tm*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-tm*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = kR*k13*(E2/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A3term = A3term1+A3term2+A3term3
  } else A3term = NULL

  if(returncpt %in% c("all", "cpt4")){
    Amterm1 = Amlast*exp(-tm*kme) +km*A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(kme-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-tm*kme)*(E2-kme)*(E3-kme)/((lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm2 = km*(exp(-tm*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-kme))+exp(-tm*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-kme))+exp(-tm*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-kme))-exp(-tm*kme)*(B*kme-C)/((lambda1-kme)*(kme-lambda2)*(kme-lambda3)))
    Amterm3 = km*kR*((E2*E3)/(lambda1*lambda2*lambda3*kme)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(kme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))-exp(-tm*kme)*(E2-kme)*(E3-kme)/(kme*(lambda1-kme)*(lambda2-kme)*(lambda3-kme)))
    Amterm = Amterm1+Amterm2+Amterm3
  } else Amterm = NULL

  return(rbind(A1term/v1, A2term/v2, A3term/v3, Amterm/v4))

}
class(pkmod3cptm) <- "pkmod"







