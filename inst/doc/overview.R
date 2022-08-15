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

## ---- libraries---------------------------------------------------------------
library(tci)
library(ggplot2)   # ggplot for plotting
library(gridExtra) # arrangeGrob to arrange plots
library(reshape2)  # melt function

## ---- echo=FALSE, eval = FALSE------------------------------------------------
#  old <- theme_set(theme_bw())
#  ggplot <- function(...) ggplot2::ggplot(...) +
#    scale_color_brewer(palette="Pastel1")
#    # scale_color_manual(values = c("black","steelblue","seagreen"))

## -----------------------------------------------------------------------------
# 1-compartment model
(mod1cpt <- pkmod(pars_pk = c(cl = 10, v = 15)))
# 3-compartment model with effect site
(mod3ecpt <- pkmod(pars_pk = c(cl = 10, q2 = 2, q3 =20, v = 15, v2 = 30, v3 = 50, ke0 = 1.2)))

## -----------------------------------------------------------------------------
# acceptable parameter names
list_parnms()

## -----------------------------------------------------------------------------
update(mod3ecpt, pars_pk = c(ke0 = 0.9), init = c(1,0.2,0.3,1))

## -----------------------------------------------------------------------------
# single infusion
(single_inf <- inf_manual(inf_tms = 0, duration = 0.5, inf_rate = 100))
# multiple infusions
(multi_inf <- inf_manual(inf_tms = c(0,3,6), duration = c(1,0.5,0.25), inf_rate = 100))

## -----------------------------------------------------------------------------
# plasma targeting for one-compartment model
inf_1cpt <- inf_tci(target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), 
                    pkmod = mod1cpt, type = "plasma")
head(inf_1cpt)

# effect-site targeting for three-compartment effect site model
inf_3ecpt <- inf_tci(target_vals = c(2,3,4,4), target_tms = c(0,2,3,10), 
                     pkmod = mod3ecpt, type = "effect")
head(inf_3ecpt)

## -----------------------------------------------------------------------------
# prediction/observation times
tms_pred <- seq(0,10,0.01)
tms_obs <- c(0.5,1,2,4,6,10)

pre <- predict(mod3ecpt, inf = inf_3ecpt, tms = tms_pred)
obs <- simulate(mod3ecpt, seed = 1, inf = inf_3ecpt, tms = tms_obs, sigma_mult = 0.2)

# plot results
dat <- data.frame(time = tms_pred, `plasma (3 cmpt)` = pre[,"c1"], 
                  `effect (ke0=1.2)` = pre[,"c4"],
                  check.names = FALSE)
datm <- melt(dat, id = "time")
dat_obs <- data.frame(time = tms_obs, con = obs, variable = "plasma (3 cmpt)")

p <- ggplot(datm, aes(x = time, y = value, color = variable)) + 
  geom_line() + 
  geom_point(data = dat_obs, aes(x = time, y = con)) +
  xlab("Minutes") + ylab("Concentration (mg/L)")
p

## -----------------------------------------------------------------------------
# evaluate with different ke0 parameter
pre_misspec <- predict(mod3ecpt, inf = inf_3ecpt, tms = tms_pred, 
                       pars_pk = c(ke0 = 0.8))
dat_misspec <- data.frame(pre_misspec, variable = "effect (ke0=0.8)", time = tms_pred)
p + geom_line(data = dat_misspec, aes(x = time, y = c4, color = variable))

## -----------------------------------------------------------------------------
# predicted concentrations
pre_1cpt <- predict(mod1cpt, inf = inf_3ecpt, tms = tms_pred)
dat_1cpt <- data.frame(pre_1cpt, variable = "plasma (1 cmpt)", time = tms_pred)
# simulated observations
obs_1cpt <- simulate(mod1cpt, seed = 1, inf = inf_3ecpt, tms = tms_obs, sigma_mult = 0.2)

p + geom_line(data = dat_1cpt, aes(x = time, y = c1, color = variable)) +
  geom_point(data = data.frame(time = tms_obs, con = obs_1cpt, variable = "plasma (1 cmpt)"), 
           aes(x = time, y = con), inherit.aes = FALSE, color = "green4")

## -----------------------------------------------------------------------------
modpd <- update(mod3ecpt, pdfn = emax, pdinv = emax_inv, 
                 pars_pd = c(e0 = 100, emx = 100, c50 = 3.5, gamma = 2.2))

## -----------------------------------------------------------------------------
inf_pd <- inf_tci(target_vals = c(70,60,50,50), target_tms = c(0,2,3,10), pkmod = modpd, type = "effect")

## -----------------------------------------------------------------------------
# predict responses
pre_pd <- predict(modpd, inf = inf_pd, tms = tms_pred)
# pd observations: 10 sec = 1/6 min
tms_pd_obs <- seq(1/6,10,1/6) 
# simulate responses with additive error and parameter misspecification
obs_pd <- simulate(modpd, seed = 1, inf = inf_pd, tms = tms_pd_obs, sigma_add = 5, 
                   pars_pk = c(ke0 = 0.7), pars_pd = c(c50 = 3, gamma = 1.8))

# plot results
dat_pd <- data.frame(time = tms_pred, `plasma (3 cmpt)` = pre_pd[,"c1"], 
                  `effect (ke0=1.2)` = pre_pd[,"c4"],
                  BIS = pre_pd[,"pdresp"],
                  check.names = FALSE)
dat_pdm <- melt(dat_pd, id = "time")
dat_pdm$type <- as.factor(ifelse(dat_pdm$variable == "BIS", "PD","PK"))
dat_pd_obs <- data.frame(time = tms_pd_obs, BIS = obs_pd, 
                         type = factor("PD"), variable = "BIS")
levels(dat_pdm$type) <- levels(dat_pd_obs$type) <- c("Bispectral Index", "Concentration (mg/L)")

ggplot(dat_pdm, aes(x = time, y = value, color = variable)) + 
  facet_wrap(type~., scales = "free", nrow = 2) +
  geom_line() + 
  geom_point(data = dat_pd_obs, aes(x = time, y = BIS)) + 
  xlab("Minutes") + ylab("")

## -----------------------------------------------------------------------------
mod_true  <- update(mod3ecpt, pars_pk = c(cl = 20, q2 = 1.5, ke0 = 1.8))
sim_ol <- simulate_tci(pkmod_prior = mod3ecpt, 
                       pkmod_true = mod_true, 
                       target_vals = c(2,3,4,4), 
                       target_tms = c(0,2,3,24),
                       obs_tms = c(1,2,3,4,8,12),
                       seed = 1)
ggplot(melt(sim_ol$resp, id.vars = c("time","type"))) + 
  geom_line(aes(x = time, y = value, color = variable)) + 
  geom_point(data = sim_ol$obs, aes(x = time, y = obs)) +
  facet_wrap(~type) +
  labs(x = "Hours", y = "Concentration (mg/L)")

## -----------------------------------------------------------------------------
mod3ecpt <- update(mod3ecpt, sigma_mult = 0.2, 
                   Omega = matrix(diag(c(1.2,0.6,1.5,0.05)), 4,4, 
                                  dimnames = list(NULL, c("cl","q2","v","ke0"))))
sim_cl <- simulate_tci(pkmod_prior = mod3ecpt, 
                       pkmod_true = mod_true, 
                       target_vals = c(2,3,4,4), 
                       target_tms = c(0,2,3,24),
                       obs_tms = c(1,2,3,4,8,12),
                       update_tms = c(6,12,16),
                       delay = 0,
                       seed = 1)
ggplot(melt(sim_cl$resp, id.vars = c("time","type"))) + 
  geom_line(aes(x = time, y = value, color = variable)) + 
  geom_point(data = sim_cl$obs, aes(x = time, y = obs)) +
  facet_wrap(~type) +
  labs(x = "Hours", y = "Concentration (mg/L)")

