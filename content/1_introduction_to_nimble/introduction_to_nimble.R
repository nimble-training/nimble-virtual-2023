## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE,
                      tidy.opts = list(width.cutoff = 60),
                      tidy = TRUE)
has_ggplot2 <- require(ggplot2)
has_mcmcplots <- require(mcmcplots)
has_coda <- require(coda)
generate_original_results <- TRUE


## -----------------------------------------------------------------------------
library(nimble)


## -----------------------------------------------------------------------------
DeerEcervi <- read.table(file.path('..', 'examples', 'DeerEcervi', 'DeerEcervi.txt'), header = TRUE)
summary(DeerEcervi)

## Create presence/absence data from counts.
DeerEcervi$Ecervi_01 <- DeerEcervi$Ecervi
DeerEcervi$Ecervi_01[DeerEcervi$Ecervi>0] <- 1
## Set up naming convention for centered and uncentered lengths for exercises later
DeerEcervi$unctrLength <- DeerEcervi$Length
## Center Length for better interpretation
DeerEcervi$ctrLength <- DeerEcervi$Length - mean(DeerEcervi$Length)
## Make a factor version of Sex for plotting
DeerEcervi$fSex <- factor(DeerEcervi$Sex)
## Make a factor and id version of Farm
DeerEcervi$fFarm <- factor(DeerEcervi$Farm)
DeerEcervi$farm_ids <- as.numeric(DeerEcervi$fFarm)


## ----eval=has_ggplot2---------------------------------------------------------
ggplot(data = DeerEcervi, 
        mapping = aes(x = ctrLength, y = Ecervi_01, color = fSex)) + 
  geom_point() + 
  geom_jitter(width = 0, height = 0.1) + 
  facet_wrap(~Farm)


## -----------------------------------------------------------------------------
DEcode <- nimbleCode({
  for(i in 1:2) {
    # Priors for intercepts and length coefficients for sex = 1 (male), 2 (female)
    sex_int[i] ~ dnorm(0, sd = 1000)
    length_coef[i] ~ dnorm(0, sd = 1000)
  }

  # Priors for farm random effects and their standard deviation.
  farm_sd ~ dunif(0, 20)
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(0, sd = farm_sd)
  }

  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    logit(disease_probability[i]) <-
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})


## -----------------------------------------------------------------------------
DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    length = DeerEcervi$ctrLength,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)

DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants)


## -----------------------------------------------------------------------------
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))
# This sets the values and *flags the nodes as data*.
DEinits <- function() {
  list(sex_int = c(0, 0),
       length_coef = c(0, 0),
       farm_sd = 1,
       farm_effect = rnorm(24, 0, 1) )
}

set.seed(123)
DEmodel$setInits(DEinits())


## -----------------------------------------------------------------------------
DEmcmc <- buildMCMC(DEmodel, enableWAIC = TRUE)


## -----------------------------------------------------------------------------
cDEmodel <- compileNimble(DEmodel) 
# First call to compileNimble in a session is slower than later calls.
cDEmcmc <- compileNimble(DEmcmc, project = DEmodel)


## -----------------------------------------------------------------------------
DEresults <- runMCMC(cDEmcmc, niter=11000, nburnin=1000, WAIC=TRUE)


## -----------------------------------------------------------------------------
# Samples
samples1 <- DEresults$samples


## -----------------------------------------------------------------------------
# WAIC (Note: there are different flavors of WAIC that can be chosen earlier.)
WAIC <- DEresults$WAIC


## ----eval=FALSE---------------------------------------------------------------
## # Run this code if you want to generate your own results.
## # They won't over-write results that come with these slides.
## library(mcmcplots)
## mcmcplot(samples1, dir = ".", filename = "Ecervi_samples_mcmcplot")


## -----------------------------------------------------------------------------
WAIC


## ----echo=FALSE, eval=(has_mcmcplots & generate_original_results)-------------
# Run the previous code to generate your own results.
library(mcmcplots)
mcmcplot(samples1, dir = ".", filename = "orig_Ecervi_samples_mcmcplot")


## ----eval = FALSE-------------------------------------------------------------
## # We haven't provided coda figures, but you can make make them if you want.
## library(coda)
## pdf("Ecervi_samples_coda.pdf")
## plot(as.mcmc(samples1))
## dev.off()


## -----------------------------------------------------------------------------
set.seed(123)
DEdataAndConstants <- c(DEconstants, 
                        list(Ecervi_01 = DeerEcervi$Ecervi_01))
samples2 <- nimbleMCMC(DEcode,
                       constants = DEdataAndConstants,
                       inits = DEinits,
                       niter = 10000,
                       nburnin = 1000,
                       nchains = 2,
                       samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
summary(samples2$samples) ## from coda
samples2$WAIC


## -----------------------------------------------------------------------------
cDEmcmc$run(niter = 11000, nburnin=1000)
samples3 <- as.matrix(cDEmcmc$mvSamples)
summary(samples3)
cDEmcmc$getWAIC()

