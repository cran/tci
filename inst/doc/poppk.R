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

## ---- libraries, echo = FALSE-------------------------------------------------
library(tci)
library(ggplot2)   # ggplot for plotting

## ---- echo=FALSE, eval = FALSE------------------------------------------------
#  old <- theme_set(theme_bw())
#  ggplot <- function(...) ggplot2::ggplot(...) +
#    scale_color_brewer(palette="Pastel1")

## -----------------------------------------------------------------------------
# create a data frame of patient covariates
data <- data.frame(ID = 1:5, AGE = seq(20,60,by=10), 
                   TBW = seq(60,80,by=5), HGT = seq(150,190,by=10), 
                   MALE = c(TRUE,TRUE,FALSE,FALSE,FALSE))
# create population PK model
pkpd_elvd <- poppkmod(data = data, drug = "ppf", model = "eleveld")

## -----------------------------------------------------------------------------
set.seed(1)
pkpd_elvd_iiv <- sample_iiv(pkpd_elvd)

set.seed(1)
pkpd_elvd_iiv2 <- poppkmod(data = data, drug = "ppf", model = "eleveld", sample = TRUE)

identical(pkpd_elvd_iiv, pkpd_elvd_iiv2)

## -----------------------------------------------------------------------------
target_vals = c(75,60,50,50)
target_tms = c(0,3,6,10)

# effect-site targeting
inf_poppk <- inf_tci(pkpd_elvd, target_vals, target_tms, "effect")
head(inf_poppk)

## -----------------------------------------------------------------------------
predict(pkpd_elvd, inf = inf_poppk, tms = c(1,3))
set.seed(1)
simulate(pkpd_elvd_iiv, nsim = 3, inf = inf_poppk, tms = c(1,3), resp_bounds = c(0,100))

## -----------------------------------------------------------------------------
# values are in terms of minutes. 1/6 = 10 seconds
obs_tms <- seq(1/6,10,1/6)

sim_ol <- simulate_tci(pkmod_prior = pkpd_elvd, pkmod_true = pkpd_elvd_iiv, 
             target_vals, target_tms, obs_tms, type = "effect", seed = 1)

## -----------------------------------------------------------------------------
plot(sim_ol)

## -----------------------------------------------------------------------------
plot(sim_ol, yvar = "c4", id = c(1,3,5), show_inf = TRUE, wrap_id = TRUE)

## -----------------------------------------------------------------------------
sim_cl <- simulate_tci(pkmod_prior = pkpd_elvd, pkmod_true = pkpd_elvd_iiv, 
             target_vals, target_tms, obs_tms, update_tms = 1:10, delay = 1/3,
               type = "effect", seed = 1)

## -----------------------------------------------------------------------------
plot(sim_cl) + 
  xlab("Minutes") + 
  ylab("Bispectral Index") + 
  ggtitle("Closed-loop simulation of Eleveld propofol model", 
          subtitle = "Minute updates, processing delay of 20 seconds")

## -----------------------------------------------------------------------------
plot(sim_cl, yvar = "c4", id = c(1,3,5), show_inf = TRUE, wrap_id = TRUE)

