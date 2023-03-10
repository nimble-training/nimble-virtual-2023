## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
has_nimbleEcology <- require(nimbleEcology)
if(!has_nimbleEcology)
  message("This module will use nimbleEcology, which you don't have installed.")
library(nimble)
library(compareMCMCs)
DO_POSTERIOR_PREDICTION <- FALSE
runCases <- TRUE
DO_PLOT <- TRUE
save_dir <- file.path("..","..","..","..","Jan2023_SwissGreatTits_saved_results")
# if needed:
# dir.create(save_dir)


## -----------------------------------------------------------------------------
if(!exists("DO_PLOT"))
    DO_PLOT <- FALSE

library(AHMbook)
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
str(SwissTits)
SwissTits$species  # Available species

# Select Great tit and covariate data from 2013 and
#   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
str(y)  # Check the matrix of count data
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

# Plot observed data: counts vs survey date (Fig. 6-9)
if(DO_PLOT) matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey data (1 = April 1)", ylab = "Count of Great Tits")

# Load unmarked, create unmarked data frame and inspect result
library(unmarked)
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
                           siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]), iLength=1/tits[,"rlength"]),
                           obsCovs=list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)                            # Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) # Summarize max counts

elev <- umf@siteCovs$elev   ;   elev2 <- elev^2
forest <- umf@siteCovs$forest   ;   forest2 <- forest^2
date <- matrix(umf@obsCovs$date, ncol = 3, byrow = TRUE)
dur <- matrix(umf@obsCovs$dur, ncol = 3, byrow = TRUE)
date[is.na(date)] <- 0   ;   date2 <- date^2
dur[is.na(dur)] <- 0   ;   dur2 <- dur^2
iRoute <- umf@siteCovs$iLength

# Design matrix for abundance model (no intercept)
lamDM <- model.matrix(~ elev + elev2 + forest + forest2 + elev:forest + elev:forest2 + iRoute)[,-1]

# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
SGT_inits <- function(){ list(N = Nst, 
                              beta0 = 0, 
                              mean.p = rep(0.5, 3), 
                              beta = runif(7, 0, 0), 
                              alpha = runif(13, 0, 0)
                              )}


# Bundle data and choose to fit simple ZIP model (model 1)
SGT_data1 <- list(y = y, 
                  lamDM = lamDM, 
                  elev = elev, 
                  date = date, 
                  dur = dur, 
                  elev2 = elev2,
                  date2 = date2, 
                  dur2 = dur2, 
                  e = 1e-06, 
                  hlam.on = 1,      # Three toggles for model components - discusss
                  hp.site.on = 1,
                  hp.survey.on = 1,
                  nsite = 263,
                  nrep = 3)

save(SGT_data1, Nst, SGT_inits, file = file.path("..", "examples", "SwissGreatTits", "prepared_data.RData"))


## -----------------------------------------------------------------------------
load(file = file.path("..", "examples", "SwissGreatTits", "prepared_data.RData"))


## -----------------------------------------------------------------------------
Section6p11_code <- nimbleCode( {
  
  # Specify priors
  # zero-inflation/suitability
  phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)
  
  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
  
  # detection
  for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1) # p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
  
  # ZIP model for abundance
  for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    N[i] ~ dpois(mu.poisson[i])
  }
  
  # Measurement error model
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # Stabilize logit
      lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
        alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
        alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
        alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
        alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
        alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
        alpha[13] * date[i,j] * dur2[i,j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
      eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
  if(DO_POSTERIOR_PREDICTION) {
    # Posterior predictive distributions of chi2 discrepancy
    for (i in 1:nsite) {
      for (j in 1:nrep) {
        y.sim[i,j] ~ dbin(p[i,j], N[i]) # Create new data set under model
        e.count[i,j] <- N[i] * p[i,j]   # Expected datum
        # Chi-square discrepancy for the actual data
        chi2.actual[i,j] <- pow((y[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
        # Chi-square discrepancy for the simulated ('perfect') data
        chi2.sim[i,j] <- pow((y.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
        # Add small value e to denominator to avoid division by zero
      }
    }
    # Add up individual chi2 values for overall fit statistic
    fit.actual <- sum(chi2.actual[1:263, 1:3])  # Fit statistic for actual data set
    fit.sim <- sum(chi2.sim[1:263, 1:3])        # Fit statistic for a fitting model
    bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
    c.hat <- fit.actual/fit.sim        # c-hat estimate
    
    # Derived parameters: Total abundance at 263 sampled sites
    Ntotal263 <- sum(N[1:263])
  }
}
)


## ----first-model--------------------------------------------------------------
DO_POSTERIOR_PREDICTION <- FALSE
m <- nimbleModel(Section6p11_code,
                 constants = SGT_data1,
                 inits = SGT_inits())
## What is the warning about?
m$initializeInfo()
m$phi
head(m$eps.p.survey)


## -----------------------------------------------------------------------------
SGT_inits_full<- function(){
    ans <- list(N = Nst, 
                beta0 = 0, 
                mean.p = rep(0.5, 3), 
                beta = runif(7, 0, 0), 
                alpha = runif(13, 0, 0)
              , phi = 0.5
              , sd.lam = 1
              , sd.p.site = 1
              , sd.p.survey = 1
              , a = rbinom(SGT_data1$nsite, prob = 0.5, size = 1)
              , eps.lam = rnorm(SGT_data1$nsite)
              , eps.p.site = rnorm(SGT_data1$nsite)
              , eps.p.survey = matrix(rnorm(SGT_data1$nsite * SGT_data1$nrep),
                                      nrow = SGT_data1$nsite)
                )
    ans$a[ ans$N > 0 ] <- 1
    ans
}


## -----------------------------------------------------------------------------
m1 <- nimbleModel(Section6p11_code,
                 constants = SGT_data1,
                 inits = SGT_inits_full())
m1$initializeInfo()
## The only nodes left not initialized are parts of y
m1$y
## The missing data are ok: they are really missing y values from rep 3.
## Note, in this case we don't need these sampled, but they will be.


## ----comp-cost-1, eval = TRUE-------------------------------------------------
MCMC1 <- buildMCMC(m1)
compiled1 <- compileNimble(m1, MCMC1)
## See that it runs for small number of iterations  
compiled1$MCMC1$run(10)
## We saw reports of an initialization problem above,
## so let's do some manual calculation to check
compiled1$m1$calculate()
## We see the compiled model is now fully
## initialized and has a valid total logProb.
m1$calculate()
## We see the uncompiled model (with original inits)
## does not give a valid logProb.
m1$logProb_y
## We see the problem is missing parts of y
m1$calculate("y")
## Indeed, calculating logProbs just of y gives NA
compiled1$m1$calculate("y")
## But again the compiled model is ok after
## MCMC initialization.

## Conclusion: What we see is that elements of y that are missing
## data make the model start with invalid log probability.
## That may not be a problem.  Those nodes will be
## initialized from their priors, so if those give decent
## values, there will be no problem.
## more later.  For now, we will move on.

## Note that when code encounters NAs or NaNs,
## it may run slower.


## ---- eval=runCases, echo=TRUE------------------------------------------------
comparison1 <- compareMCMCs(modelInfo = list(model = m1),
                            MCMCcontrol = list(niter = 120000, nburnin = 20000, thin = 100))
save(comparison1, file = file.path(save_dir, "ZIPNmix_case1.Rdata"))


## ---- eval=!runCases, echo=FALSE----------------------------------------------
## load(file.path(save_dir, "ZIPNmix_case1.Rdata"))


## -----------------------------------------------------------------------------
# How long did 120000 iterations take?
comparison1$nimble$times$sampling
# Mixing and efficiency
comparison1$nimble$metrics


## ---- echo=FALSE--------------------------------------------------------------
Section6p11_code_jags <- nimbleCode( {
  
  # Specify priors
  # zero-inflation/suitability
  phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)
  
  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
  
  # detection
  for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1) # p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
  
  # ZIP model for abundance
  for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    N[i] ~ dpois(mu.poisson[i])
  }
  
  # Measurement error model
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # Stabilize logit
      lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
        alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
        alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
        alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
        alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
        alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
        alpha[13] * date[i,j] * dur2[i,j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
      eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
})


## ---- eval=runCases, echo=TRUE------------------------------------------------
comparison1jags <- compareMCMCs(modelInfo =
                                  list(code = Section6p11_code_jags, constants = SGT_data1, inits = SGT_inits_full()),
                                MCMCcontrol = list(niter = 120000, nburnin = 20000, thin = 100),
                                MCMCs = "jags")
save(comparison1jags, file = file.path(save_dir, "ZIPNmix_case1jags.Rdata"))


## ---- eval=!runCases, echo=FALSE----------------------------------------------
## load(file.path(save_dir, "ZIPNmix_case1jags.Rdata"))


## -----------------------------------------------------------------------------
# How long did 120000 iterations take?
comparison1jags$jags$times$sampling
# Mixing and efficiency
comparison1jags$jags$metrics


## -----------------------------------------------------------------------------
Section6p11_code_grouped <- nimbleCode( {
  # Specify priors
  # zero-inflation/suitability
  phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)
  
  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
  
  # detection
  for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1) # p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
  
  # ZIP model for abundance
  # new
  loglam.fixed[1:nsite] <- beta0 + (lamDM[1:nsite, 1:7] %*% beta[1:7])[1:nsite,1]
  for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    loglam[i] <- loglam.fixed[i] + eps.lam[i] * hlam.on
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    N[i] ~ dpois(mu.poisson[i])
  }
  # Measurement error model
  # new
  for(j in 1:nrep) {
    lp.fixed[1:nsite,j] <- alpha0[j] + alpha[1] * elev[1:nsite] + alpha[2] * elev2[1:nsite] +
      alpha[3] * date[1:nsite,j] + alpha[4] * date2[1:nsite,j] +
      alpha[5] * dur[1:nsite,j] + alpha[6] * dur2[1:nsite,j] +
      alpha[7] * elev[1:nsite] * date[1:nsite,j] + alpha[8] * elev2[1:nsite] * date[1:nsite,j] +
      alpha[9] * elev[1:nsite] * dur[1:nsite,j] + alpha[10] * elev[1:nsite] * dur2[1:nsite,j] +
      alpha[11] * elev2[1:nsite] * dur[1:nsite,j] + alpha[12] * date[1:nsite,j] * dur[1:nsite,j] +
      alpha[13] * date[1:nsite,j] * dur2[1:nsite,j]
  }
  
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # Stabilize logit

      ## new
      lp[i,j] <- lp.fixed[i, j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
      
      eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
}
)


## ----comp-cost-2, eval=runCases-----------------------------------------------
m2 <- nimbleModel(Section6p11_code_grouped,
                 constants = SGT_data1,
                 inits = SGT_inits_full())
comparison2 <- compareMCMCs(modelInfo = list(model = m2),
                            MCMCcontrol = list(niter = 120000, nburnin = 20000, thin = 100))
save(comparison2, file = file.path(save_dir, "ZIPNmix_case2.Rdata"))


## ---- eval=!runCases, echo=FALSE----------------------------------------------
## load(file.path(save_dir, "ZIPNmix_case2.Rdata"))


## -----------------------------------------------------------------------------
# How long did 120000 iterations take?
comparison2$nimble$times$sampling
# Mixing and efficiency
comparison2$nimble$metrics


## ---- eval=FALSE--------------------------------------------------------------
## y[i, 1:num_surveys] ~ dNmixture_v(lambda[i],
##                                  prob[1:num_surveys],
##                                  Nmin = -1, Nmax = -1,
##                                  len = num_surveys)


## -----------------------------------------------------------------------------
Section6p11_code_dNmixture <- nimbleCode( {
  phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)
  
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
  
  # detection
  for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1) # p intercept for occasions 1-3
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
  
  # ZIP model for abundance
  # new
  loglam.fixed[1:nsite] <- beta0 + (lamDM[1:nsite, 1:7] %*% beta[1:7])[1:nsite,1]
  for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    
    ## original
    ## loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
    ## new
    loglam[i] <- loglam.fixed[i] + eps.lam[i] * hlam.on
    
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    # N[i] ~ dpois(mu.poisson[i])
  }
  
  # Measurement error model
  # new
  for(j in 1:nrep) {
    lp.fixed[1:nsite,j] <- alpha0[j] + alpha[1] * elev[1:nsite] + alpha[2] * elev2[1:nsite] +
      alpha[3] * date[1:nsite,j] + alpha[4] * date2[1:nsite,j] +
      alpha[5] * dur[1:nsite,j] + alpha[6] * dur2[1:nsite,j] +
      alpha[7] * elev[1:nsite] * date[1:nsite,j] + alpha[8] * elev2[1:nsite] * date[1:nsite,j] +
      alpha[9] * elev[1:nsite] * dur[1:nsite,j] + alpha[10] * elev[1:nsite] * dur2[1:nsite,j] +
      alpha[11] * elev2[1:nsite] * dur[1:nsite,j] + alpha[12] * date[1:nsite,j] * dur[1:nsite,j] +
      alpha[13] * date[1:nsite,j] * dur2[1:nsite,j]
  }
  
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    y[i, 1:nrep_by_site[i]] ~ dNmixture_v(mu.poisson[i], p[i, 1:nrep_by_site[i]], Nmin = -1, Nmax = -1, len = nrep_by_site[i])
    p[i, 1:nrep_by_site[i]] <- 1 / (1 + exp(-lp.lim[i,1:nrep_by_site[i]]))
    for (j in 1:nrep_by_site[i]){
      #        y[i,j] ~ dbin(p[i,j], N[i])
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # Stabilize logit
      ## new
      lp[i,j] <- lp.fixed[i, j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
      
      eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
})


## ----comp-cost-4, eval=runCases-----------------------------------------------
nrep_by_site <- apply(y, 1, function(x) sum(!is.na(x)))
SGT_data1$nrep_by_site <- nrep_by_site
m3 <- nimbleModel(Section6p11_code_dNmixture,
                 constants = SGT_data1,
                 inits = SGT_inits_full())
comparison3 <- compareMCMCs(modelInfo = list(model = m3),
                            MCMCcontrol = list(niter = 120000, nburnin = 20000, thin = 100))
save(comparison3, file = file.path(save_dir, "ZIPNmix_case3.Rdata"))


## ---- eval=!runCases, echo=FALSE----------------------------------------------
## load(file.path(save_dir, "ZIPNmix_case3.Rdata"))


## -----------------------------------------------------------------------------
# How long did 120000 iterations take?
comparison3$nimble$times$sampling
# Mixing and efficiency
comparison3$nimble$metrics


## -----------------------------------------------------------------------------
Section6p11_code_grouped_dlogis <- nimbleCode( {
  
  # Specify priors
  # zero-inflation/suitability
  phi ~ dunif(0,1)          # proportion of suitable sites (probability of being not a structural 0)
  theta <- 1-phi            # zero-inflation (proportion of unsuitable)
  ltheta <- logit(theta)
  
  # abundance
  beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
  for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
  
  # detection
  for(j in 1:3){
    alpha0[j] ~ dlogis(0,1)
  }
  for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
  }
  tau.p.site <- pow(sd.p.site, -2)
  sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
  tau.p.survey <- pow(sd.p.survey, -2)
  sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
  
  # ZIP model for abundance
  # new
  loglam.fixed[1:nsite] <- beta0 + (lamDM[1:nsite, 1:7] %*% beta[1:7])[1:nsite,1]
  for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    
    ## original
    ## loglam[i] <- beta0 + inprod(beta[1:7], lamDM[i, 1:7]) + eps.lam[i] * hlam.on
    ## new
    loglam[i] <- loglam.fixed[i] + eps.lam[i] * hlam.on
    
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    N[i] ~ dpois(mu.poisson[i])
  }
  
  # Measurement error model
  # new
  for(j in 1:nrep) {
    lp.fixed[1:nsite,j] <- alpha0[j] + alpha[1] * elev[1:nsite] + alpha[2] * elev2[1:nsite] +
      alpha[3] * date[1:nsite,j] + alpha[4] * date2[1:nsite,j] +
      alpha[5] * dur[1:nsite,j] + alpha[6] * dur2[1:nsite,j] +
      alpha[7] * elev[1:nsite] * date[1:nsite,j] + alpha[8] * elev2[1:nsite] * date[1:nsite,j] +
      alpha[9] * elev[1:nsite] * dur[1:nsite,j] + alpha[10] * elev[1:nsite] * dur2[1:nsite,j] +
      alpha[11] * elev2[1:nsite] * dur[1:nsite,j] + alpha[12] * date[1:nsite,j] * dur[1:nsite,j] +
      alpha[13] * date[1:nsite,j] * dur2[1:nsite,j]
  }
  
  for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # Stabilize logit
      ## new
      lp[i,j] <- lp.fixed[i, j] +
        eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
      
      eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
  }
}
)



## ---- eval=runCases-----------------------------------------------------------
assignRWB <- function(model) {
  MCMCconf <- configureMCMC(model)
  MCMCconf$removeSamplers(c("beta0","beta"))
  # Set up block for one set of coefficients, with repetition
  MCMCconf$addSampler(c("beta0","beta"), type = "RW_block", tries = 5, adaptInterval=20)
  MCMCconf$removeSamplers(c("alpha0", "alpha"))
  # Set up block for the other set of coefficients, with repetition
  MCMCconf$addSampler(c("alpha0", "alpha"), type = "RW_block", tries = 5, adaptInterval=20)
  MCMCconf
}

m4 <- nimbleModel(Section6p11_code_grouped_dlogis,
                 constants = SGT_data1,
                 inits = SGT_inits_full())

comparison4 <- compareMCMCs(modelInfo = list(model = m4),
                            MCMCcontrol = list(niter = 120000, nburnin = 20000, thin = 100),
                            nimbleMCMCdefs = list(RWB = assignRWB),
                            MCMCs = "RWB",
                            conversions =
                              list(RWB = list(
                                "mean.p[1]" = "expit(`alpha0[1]`)",
                                "alpha0[1]" = NULL,
                                "mean.p[2]" = "expit(`alpha0[2]`)",
                                "alpha0[2]" = NULL,
                                "mean.p[3]" = "expit(`alpha0[3]`)",
                                "alpha0[3]" = NULL)))
save(comparison4, file = file.path(save_dir, "ZIPNmix_case4.Rdata"))


## ---- eval=!runCases, echo=FALSE----------------------------------------------
## load(file.path(save_dir, "ZIPNmix_case4.Rdata"))


## -----------------------------------------------------------------------------
# How long did 120000 iterations take?
comparison4$RWB$times$sampling
# Mixing and efficiency
comparison4$RWB$metrics


## -----------------------------------------------------------------------------
c1 <- renameMCMC(comparison1, "original", "nimble")
c2 <- renameMCMC(comparison2, "grouped", "nimble")
c3 <- renameMCMC(comparison3, "Nmix", "nimble")
make_MCMC_comparison_pages(c(comparison1jags, c1, c2, c3, comparison4),
                           dir = "Swiss_Great_Tit_MCMC_comparisons")

