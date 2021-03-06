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

## ---- load-data---------------------------------------------------------------
library(tci)

data("eleveld_pk")
data("eleveld_pd")

# 122 patients had both pk and pd measurements
pkpd_eb <- merge(eleveld_pk, eleveld_pd)

## ---- plc-sim-1---------------------------------------------------------------
# Evaluate Eleveld model at covariate values for patient 1
prior_pars_id1 <- eleveld_poppk(df = pkpd_eb[1,]) 
prior_pars_id1

# tci infusions are calculated using prior parameter estimates
pars_pk_id1 <- unlist(prior_pars_id1[,c("V1","V2","V3","CL","Q2","Q3","KE0")])
# BIS0 is repeated since BIS at maximum effect = BIS0 for Eleveld model
pars_pd_id1 <- unlist(prior_pars_id1[,c("CE50","GAMMA","GAMMA2","BIS0","BIS0")])

## ---- olc-sim-2---------------------------------------------------------------
olc_inf_id1 <- tci_pd(pdresp = c(50,50), # target BIS = 50
                      tms = c(0,5), # target for 5 minutes
                      pkmod = pkmod3cptm, # pk model
                      pdmod = emax_eleveld, # pd model
                      pdinv = inv_emax_eleveld, # inverse pd model
                      pars_pk = pars_pk_id1,
                      pars_pd = pars_pd_id1)
plot(olc_inf_id1)

## ---- sim-data-olc------------------------------------------------------------
# true set of patient parameters given by EB estimates
eb_pars_pk_id1 <- unlist(pkpd_eb[1,c("CL","Q2","Q3","V1","V2","V3","KE0")])
eb_pars_pd_id1 <- unlist(pkpd_eb[1,c("E50","GAM","GAM1","EMAX","EMAX")])

# residual error term
eb_sigma_id1 <- unlist(pkpd_eb[1,"RESD"])

# BIS delay - convert from minutes to seconds
eb_bd_id1 <- unlist(pkpd_eb[1,"ALAG1"]) * 60

# simulate response to infusions
set.seed(1)
datasim_id1 <- gen_data(inf = olc_inf_id1, 
                        pkmod = pkmod3cptm, 
                        pdmod = emax_eleveld, 
                        pars_pk0 = eb_pars_pk_id1, 
                        pars_pd0 = eb_pars_pd_id1, 
                        sigma_add = eb_sigma_id1, 
                        delay = eb_bd_id1)

# plot results
plot(datasim_id1)

## ---- closed-loop, cache=TRUE-------------------------------------------------
cl_targets <- function(time, target){
  data.frame(time = time, target = target)
}

cl_updates <- function(time, full_data = TRUE, plot_progress = FALSE){
  data.frame(time = time, full_data = full_data, plot_progress = plot_progress)
}

# 1) Data frame of targets
targets = cl_targets(time = seq(0,10,1/6),
                     target = 50)

# 2) Update times
updates <- cl_updates(time = seq(2,10,2))

# 3) prior parameter values based on point estimates at covariate values 
prior_par_list <- list(
  pars_pkpd = prior_pars_id1[c("CL","Q2","Q3","V1","V2","V3",
                                 "KE0","CE50","GAMMA","GAMMA2",
                                 "BIS0","BIS0")],
              pk_ix = 1:7, # indices of PK parameters
              pd_ix = 8:12,
              fixed_ix = 9:12, # indices of parameters that aren't updated
              err = prior_pars_id1[1,"SIGMA"], # residual error
              sig = eleveld_vcov(prior_pars_id1, 
                                 rates = FALSE)[[1]]
              )

# 4) true parameter values based on EB values
true_par_list <- list(pars_pkpd = pkpd_eb[1,c("CL","Q2","Q3","V1",
                                              "V2","V3","KE0",
                                              "E50","GAM","GAM1",
                                              "EMAX","EMAX")],
                      pk_ix = 1:7,
                      pd_ix = 8:12,
                      fixed_ix = 9:12,
                      err = pkpd_eb[1,"RESD"],
                      delay = pkpd_eb[1,"ALAG1"]) 

## -----------------------------------------------------------------------------
# run simulation
bayes_sim <- bayes_control(targets = targets, 
                          updates = updates,
                          prior = prior_par_list, 
                          true_pars = true_par_list)
plot(bayes_sim)

