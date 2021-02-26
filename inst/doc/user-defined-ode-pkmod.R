## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- mrgsolve-example, warning=FALSE-----------------------------------------
library(tci)
library(mrgsolve)

# define two-compartment model based on ODEs
code <- '
$PARAM CL=10, V=10, Q = 10, V2=10
$CMT A1 A2
$ODE
dxdt_A1 = -(CL/V)*A1 - (Q/V)*A1 + (Q/V2)*A2;
dxdt_A2 = -(Q/V2)*A2 + (Q/V)*A1;
'
mod2cpt <- mcode("example", code)

## -----------------------------------------------------------------------------
# wrapper function to evaluate pk2 in format used by tci
pkmod2cpt_ode <- function(tm, kR, pars, 
                          pkm = pk2, 
                          init = c(0,0), 
                          inittm = 0){

  # begin times at 0 and end at last time evaluated
  tm <- tm - inittm
  end_inf <- max(tm)
  
  # name parameter/initial concentrations
  names(pars) <- toupper(names(pars))
  names(init) <- c("A1","A2")
  
  # pass parameters as list
  pars <- sapply(pars, as.list)
  vols <- unlist(pars[c("V","V2")])
  
  # update parameters and initial values (as amounts)
  mod2cpt <- update(mod2cpt, param = pars, init = init*vols)
  
  # dosing regimen - mrgsolve function in terms of amount infused
  event <- ev(amt =  kR*end_inf, time = 0, tinf = end_inf)
  
  # simulate responses (skip tm=0 unless specified)
  dat <- mrgsim_q(x = mod2cpt, # pk model
                  data = event, # dosing event
                  stime = tm) # evaluation times
  
  # skip tm=0 unless specified in tm
  dat <- dat@data[-1,]
  
  # return concentrations with compartments in rows and times in columns
  cons <- t(dat[,c("A1","A2")]) / vols
  
  rownames(cons) <- colnames(cons) <- NULL
  return(cons)
}
class(pkmod2cpt_ode) <- "pkmod"

## ---- fig.height= 6, fig.width = 6--------------------------------------------
# create dose
dose <- create_intvl(
  as.matrix(cbind(time = c(0.5,4,4.5,10), 
                  infrt = c(100,0,100,0)))
)

# model parameters
pars_2cpt <- c(CL = 15, V = 10, Q = 10, V2 = 20)

# predict concentrations
predict(pkmod2cpt_ode, 
        pars = pars_2cpt, 
        inf = dose, 
        tms = c(1,2,3))

# plot concentrations
plot(pkmod2cpt_ode, 
     pars = pars_2cpt, 
     inf = dose, 
     title = "Concentrations for a user-defined 2-compartment model")

## ---- user-defined-tci, fig.height= 6, fig.width = 6--------------------------
# target concentration = 1 for times 0 - 5
# target concentration = 2 for times 5 - 10
tci_times <- c(0,5,10)
tci_targets <- c(1,2,2)

# plasma-targeting
inf_3cpt_plasma <- tci(Ct = tci_targets, 
                       tms = tci_times, 
                       pkmod = pkmod2cpt_ode, 
                       pars = pars_2cpt, 
                       tci_alg = "plasma")
# print infusion schedule
head(inf_3cpt_plasma)

# plot infusion schedule
plot(inf_3cpt_plasma)

