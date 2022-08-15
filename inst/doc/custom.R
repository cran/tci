## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, 
  message = FALSE, 
  warning = FALSE, 
  fig.align="center",
  fig.height= 4, 
  fig.width = 6
)

## ----mrgsolve-model, eval = TRUE, echo = TRUE---------------------------------
library(tci)
library(mrgsolve) # implement ODE equation
library(xtable)   # printing tables
library(ggplot2)  # plotting results
library(reshape2) # melt function

form <- '
$PARAM V1 = 7.88, V2=23.9, V3=13.8, CL1=5, CL2=0.828, CL3=0.0784,
k10 = 0.172, k12=0.373, k21=0.103, k13=0.0367, k31=0.0124
$CMT A1 A2 A3
$ODE
dxdt_A1 = k21*A2 + k31*A3 - (k12+k13+k10)*A1 - CL1/V1*A1;
dxdt_A2 = k12*A1 - k21*A2 - CL2/V2*A2;
dxdt_A3 = k13*A1 - k31*A3 - CL3/V3*A3;
'
mrg_mod_remif <- mcode("remifentanil", form)

## ----remif-wrapper------------------------------------------------------------
pk_remif <- function(tm, kR, pars, init = c(0,0,0)){

  # allow lowercase names
  names(pars) <- toupper(names(pars))
  # store volume
  vols <- pars[c("V1","V2","V3")]
  A0 <- init*vols # initial amounts
  names(A0) <- c("A1","A2","A3") # names required by mrgsolve
  
  # pass parameters as list
  pars <- sapply(pars, as.list)

  # update parameters and initial values (as amounts)
  mrg_mod_remif <- update(mrg_mod_remif, param = pars, init = A0)

  # dosing regimen - mrgsolve function in terms of amount infused
  event <- ev(amt =  kR*max(tm), time = 0, tinf = max(tm))

  # simulate responses (skip tm=0 unless specified)
  dat <- mrgsim_q(x = mrg_mod_remif, # pk model
                  data = event, # dosing event
                  stime = tm) # evaluation times

  # skip tm=0 unless specified in tm
  dat <- dat@data[-1,]

  # return concentrations with compartments in rows and times in columns
  cons <- t(dat[,c("A1","A2","A3")]) / vols

  rownames(cons) <- colnames(cons) <- NULL
  return(cons)
}

## ----Cascone-table, results='asis', echo = FALSE------------------------------
tab_cascone <- data.frame(Parameter1 = c("$V_1$","$V_2$","$V_3$","$CL_1$",
                                         "$CL_2$","$CL_3$"),
                          Optimized_val1 = c("7.88 L","23.9 L","13.8 L",
                                             "2.08 L/min","0.828 L/min","0.0784 L/min"),
                          Parameter2 = c("$k_{10}$","$k_{12}$","$k_{21}$","$k_{13}$","$k_{31}$",""),
                          Optimized_val2 = c("0.172/min","0.373/min","0.103/min","0.0367/min","0.0124/min",""))

names(tab_cascone) <- c("Parameter", "Optimized value", "Parameter", "Optimized value")
tab_cascone_df <- xtable(tab_cascone,
                         caption = "Reproduction of Table 1 from Cascone et al. (2013): Values and dimensions of the three-compartmental model parameters.",
                         label = "tab:Cascone-tab1")
print(tab_cascone_df,
      include.rownames = FALSE,
      caption.placement = "top",
      timestamp = NULL,
      booktabs = TRUE,
      hline.after = c(-1,0,6),
      sanitize.text.function = function(x){x},
      include.colnames=TRUE)

## -----------------------------------------------------------------------------
dose_remi <- inf_manual(inf_tms = 0, inf_rate = 60, duration = 20)
pars_remif <- c(V1 = 7.88, V2=23.9, V3=13.8, CL1=5, CL2=0.828, CL3=0.0784,
                k10 = 0.172, k12=0.373, k21=0.103, k13=0.0367, k31=0.0124)
mod_remif <- pkmod(pkfn = pk_remif, pars_pk = pars_remif)
p1 <- predict(mod_remif, inf = dose_remi, tms = 0:80)
ggplot(melt(data.frame(time = 0:80, p1), id = "time"), 
       aes(x = time, y = value, color = variable)) +
  geom_line()

## ----alternate-effect-site-alg, echo = TRUE-----------------------------------
tci_plasma_lim <- function(Ct, pkmod, dtm = 1/6, maxrt = 1200,
                          lim_amt = 0.5, ecmpt = NULL, tmax_search = 20, 
                          cetol = 0.05, cptol = 0.1, ...){

  pkmod <- update(pkmod,...)
  
  # if effect-site concentration is close to target,
  # switch to plasma targeting
  if(with(pkmod,(Ct - init[ecmpt]) / Ct <  cetol &
     (Ct - init[pcmpt])/Ct <= cptol))
    return(tci_plasma(Ct, pkmod = pkmod, dtm = dtm, maxrt = maxrt))
  
  # maximum tolerable plasma concentration
  Cp_max <- Ct + lim_amt
  
  # infusion required to reach Cp_max
  pinf <- tci_plasma(Ct = Cp_max, pkmod = pkmod, dtm = dtm, maxrt = maxrt)
  
  # Administer dtm-minute infusion
  unit_inf <- inf_manual(inf_tms = 0, inf_rate = pinf, duration = dtm)
  
  # Calculate maximum effect-site concentration
  CeP <- function(tm) predict(pkmod, inf = unit_inf, tms = tm)[,pkmod$ecmpt]
  Ce_max <- optimize(CeP, c(0,20), maximum = TRUE)$objective
  
  # if max Ce < Ct administer infusion to reach maximum target
  if(Ce_max <= Ct + cetol*Ct) 
    infrt <- pinf
  else
    infrt <- tci_effect_only(Ct, pkmod, dtm, maxrt = maxrt)

  return(infrt)
}

## ----base-custom-alg, echo = TRUE---------------------------------------------
mod3ecpt <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2))
tci_plasma_lim(Ct = 2, pkmod = mod3ecpt, lim_amt = 0.25)

## ----tci-custom-alg, echo = TRUE, fig.cap= "Evaluation of user-defined TCI effect-site algorithm. \\label{fig:custom-tci-alg}"----
# tci target concentrations
tci_targets <- cbind(value = c(1,2,2.5,2.5), time = c(0,3,7,10))
# calculate infusion schedule using plasma-limiting algorithm
plim_inf <- inf_tci(target_vals = c(1,2,2.5,2.5),
                    target_tms = c(0,3,7,10),
                    pkmod = mod3ecpt, 
                    custom_alg = tci_plasma_lim, 
                    lim_amt = 0.25)
head(plim_inf)

## -----------------------------------------------------------------------------
# effect-site targeting
eff_inf <- inf_tci(target_vals = c(1,2,2.5,2.5),
                    target_tms = c(0,3,7,10), 
                   pkmod = mod3ecpt, 
                   type = "effect")

## -----------------------------------------------------------------------------
# predict responses
tms_pred <- seq(0,10,0.1)
plim_pred <- predict(mod3ecpt, plim_inf, tms_pred)
eff_pred <- predict(mod3ecpt, eff_inf, tms_pred)

# plot results
dat <- data.frame(time = tms_pred, 
                  `plasma (custom)` = plim_pred[,"c1"], 
                  `effect (custom)` = plim_pred[,"c4"],
                  `plasma (effect)` = eff_pred[,"c1"], 
                  `effect (effect)` = eff_pred[,"c4"],
                  check.names = FALSE)
datm <- melt(dat, id = "time")
datm$algorithm <- ifelse(datm$variable %in% c("plasma (custom)","effect (custom)"),
                         "Plasma-limiting", "Effect-site")
ggplot(datm, aes(x = time, y = value, color = variable, linetype = algorithm)) + 
  geom_line() + 
  xlab("Minutes") + 
  ylab("Concentration (mg/L)") +
  ggtitle(label = "Plasma-limiting effect-site TCI algorithm")

