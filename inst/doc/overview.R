## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, 
  message = FALSE, 
  warning = FALSE, 
  fig.align="center",
  fig.height= 6, 
  fig.width = 6
)

## ---- libraries---------------------------------------------------------------
library(tci)
library(knitr)
library(gridExtra)

## ---- dose-object-------------------------------------------------------------
# e.g. infusion rates of 100 ug/min for 30 sec intervals at 0,4 minutes
dose <- create_intvl(
  as.matrix(cbind(time = c(0.5,4,4.5,10), 
                  infrt = c(100,0,100,0))),
  inittm = 0
)
dose

## ---- predict-mod-------------------------------------------------------------
# model parameters 
pars_3cpt <- c(k10=1.5,k12=0.15,k21=0.09,k13=0.8,
               k31=0.8,v1=10,v2=15,v3=100,ke0=1)

# predict concentrations of a three-compartment model with effect-site at
# times 1, 2, 8 minutes
predict_pkmod(pkmod3cptm, 
              pars = pars_3cpt,
              inf = dose, 
              tms = c(1,2,8))

## -----------------------------------------------------------------------------
# plot concentrations
plot(pkmod3cptm, inf = dose, pars = pars_3cpt,
     title = "Concentrations for a 3 compartment model with an effect site")

## ---- plasma-tci-1cpt, fig.height=6-------------------------------------------
# calculate infusion rates
kR_plasma <- tci_plasma(Cpt = 2, # target plasma concentration
                         dt = 1, # duration of infusion
                         pkmod = pkmod3cptm, # pk model
                         pars = pars_3cpt) # pk model parameters
kR_plasma

kR_effect <- tci_effect(Cet = 2, dt = 1, pkmod = pkmod3cptm, pars = pars_3cpt)

# set up dosing schedules
inf_2min_plasma <- create_intvl(
  data.frame(time = c(1, 20),
             infrt = c(kR_plasma,0))
  )

inf_2min_effect <- create_intvl(
  data.frame(time = c(1, 20),
             infrt = c(kR_effect,0))
  )

# plot results
plt_plasma <- plot(pkmod3cptm, 
                   pars = pars_3cpt, 
                   inf =  inf_2min_plasma, 
                   title = "Plasma targeting a concentration of 2 mg/L at 1 minute")

plt_effect <- plot(pkmod3cptm, 
                   pars = pars_3cpt, 
                   inf =  inf_2min_effect, 
                   title = "Effect-site targeting a concentration of 2 mg/L at 1 minute")

grid.arrange(plt_plasma, plt_effect)

## ---- tci-algs----------------------------------------------------------------
# target concentrations of 2 for 0-5 minutes and 3 for 5-10 minutes.
tci_times <- c(0,5,10)
tci_targets <- c(2,3,3)

# infusions for effect-site targeting
inf_3cpt_effect <- tci(Ct = tci_targets, 
                       tms = tci_times, 
                       pkmod = pkmod3cptm, 
                       pars = pars_3cpt)

plot(inf_3cpt_effect, title = "Effect-site targeting for three-compartment model")

## ---- poppk-------------------------------------------------------------------
patient_dat <- data.frame(AGE  = c(20,40,65),
                          TBM  = c(50,70,90),
                          HGT  = c(150,170,200),
                          MALE = c(TRUE,FALSE,TRUE))

# evaluate at covariate values and return clearance parameters
patient_pk <- schnider_poppk(patient_dat, rand = FALSE, rate = FALSE)
patient_pk

# evaluate TCI for patient 1
tci_patient1 <- tci(Ct = tci_targets, 
                    tms = tci_times, 
                    pkmod = pkmod3cptm, 
                    pars = patient_pk[1,], 
                    tci_alg = "effect")
head(round(tci_patient1,3))

