## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
generate_original_results <- FALSE
library(nimble)


## ---- eval=FALSE--------------------------------------------------------------
## nimbleCode({
##   # other code not shown
##   predicted[i] <- my_function(params[i, 1:p])
##   # other code not shown
## })


## ---- eval=FALSE--------------------------------------------------------------
## nimbleCode({
##   # other code not shown
##   y[i] ~ dmy_distribution(param[i], theta[i])
##   # other code not shown
## })


## ----load_dipper--------------------------------------------------------------
dipper_example_dir <- file.path("..", "..", "content", "examples","dipper")
dipper <- read.csv(file.path(dipper_example_dir,"dipper.csv"))
y <- as.matrix(dipper[ , 1:7])
first <- apply(y, 1, function(x) min(which(x !=0))) # first capture occasion
y <- y[ first != 7, ] # remove records with first capture on last occasion
head(y)


## ----dipper_basic-------------------------------------------------------------
dipper_code_basic <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  # likelihood
  for (i in 1:N){
    z[i,first[i]] <- 1
    for (t in (first[i]+1):T){
      z[i,t] ~ dbern(phi * z[i,t-1]) # z = 0 for dead, 1 for alive
      y[i,t] ~ dbern(p * z[i,t])     # y = 0 for not-captured, 1 for captured
    }}
  })


## ----setupInputs--------------------------------------------------------------
dipper_constants <- list(N = nrow(y),      # Individuals
                         T = ncol(y), # Occasions
                         first = first)    # vector of first detection occasions
zinits <- y                  # 0s and 1s
zdata <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
# Set z values for non-detections before last detection
#  to 1, so they will be treated as data.
for(i in 1:nrow(zinits)) {
  known_alive <- range(which(zinits[i,] == 1))
  zinits[i, known_alive[1] : known_alive[2] ] <- NA
  zdata[i, known_alive[1] : known_alive[2] ] <- 1
}
dipper_data <- list(y = y, z = zdata)   # 0s and 1s
dipper_inits <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
head(dipper_data$z)
head(dipper_inits()$z)


## -----------------------------------------------------------------------------
samples <- nimbleMCMC(dipper_code_basic, dipper_constants, dipper_data,
                      dipper_inits(), samplesAsCodaMCMC = TRUE)
summary(samples)


## -----------------------------------------------------------------------------
dipper_code_dCJS <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  for (i in 1:N){
    y[i, first[i]:T] ~ dCJS(phi, p)
  }
})


## -----------------------------------------------------------------------------
dCJS_R <- function (x, probSurvive, probCapture, log = FALSE) {
  probAliveGivenHistory <- 1
  logProbData <- 0
  for (t in 2:length(x)) {
    probAlive <- probAliveGivenHistory * probSurvive
    if (x[t] == 1) {
      probThisObs <- probAlive * probCapture
      probAliveGivenHistory <- 1
    } else {
      probAliveNotSeen <- probAlive * (1 - probCapture)
      probThisObs <- probAliveNotSeen + (1 - probAlive)
      probAliveGivenHistory <- probAliveNotSeen/probThisObs
    }
    logProbData <- logProbData + log(probThisObs)
  }
  if (log) return(logProbData)
  return(exp(logProbData))
}


## -----------------------------------------------------------------------------
y[5,] # A good example capture history
dCJS_R(y[5,], probSurvive = 0.7, probCapture = 0.5, log = TRUE)


## -----------------------------------------------------------------------------
dCJS <- nimbleFunction(
  run = function(x = double(1),            # vector
                 probSurvive = double(0),  # scalar
                 probCapture = double(0),  # scalar
                 log = integer(0, default = 0)) {  # integer scalar
    returnType(double())  # scalar return, can be anywhere
    probAliveGivenHistory <- 1
    logProbData <- 0
    for (t in 2:length(x)) {
      probAlive <- probAliveGivenHistory * probSurvive
      if (x[t] == 1) {
        probThisObs <- probAlive * probCapture
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture)
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen/probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) return(logProbData)
    return(exp(logProbData))
  }
)


## -----------------------------------------------------------------------------
dCJS(y[5,], probSurvive = 0.7, probCapture = 0.5, log = TRUE)


## ---- eval=FALSE--------------------------------------------------------------
## debugonce(dCJS)
## dCJS(y[5,], probSurvive = 0.7, probCapture = 0.5, log = TRUE)


## -----------------------------------------------------------------------------
dipper_code_dCJS <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  for (i in 1:N){
    y[i, first[i]:T] ~ dCJS(phi, p)
  }
})


## -----------------------------------------------------------------------------
dipper_model <- nimbleModel(code = dipper_code_dCJS,
                            constants = dipper_constants,
                            data = dipper_data,     # data can be set later.
                            inits = dipper_inits()  # inits can be set later.
                            )                       # dimensions is also a useful argument.


## -----------------------------------------------------------------------------
dipper_model$calculate()
dipper_model$calculate("y[5,]")


## ---- eval=FALSE--------------------------------------------------------------
## debugonce(dCJS)
## dipper_model$calculate("y[5,]")


## -----------------------------------------------------------------------------
dipper_MCMC <- buildMCMC(dipper_model)


## ---- eval=FALSE--------------------------------------------------------------
## debug(dCJS)
## dipper_MCMC$run(niter = 5)
## undebug(dCJS)


## -----------------------------------------------------------------------------
C_dCJS <- compileNimble(dCJS)


## -----------------------------------------------------------------------------
C_dCJS(y[5,], 0.7, 0.5, log=TRUE)


## ---- echo=FALSE--------------------------------------------------------------
dipper_model <- nimbleModel(code = dipper_code_dCJS,
                            constants = dipper_constants,
                            data = dipper_data,     # data can be set later.
                            inits = dipper_inits()  # inits can be set later.
                            )                       # dimensions is also a useful argument.


## ---- echo=FALSE--------------------------------------------------------------
dCJS <- nimbleFunction(
  run = function(x = double(1),            # vector
                 probSurvive = double(0),  # scalar
                 probCapture = double(0),  # scalar
                 log = integer(0, default = 0)) {  # integer scalar
    returnType(double())  # scalar return, can be anywhere
    probAliveGivenHistory <- 1
    logProbData <- 0
    for (t in 2:length(x)) {
      probAlive <- probAliveGivenHistory * probSurvive
      if (x[t] == 1) {
        probThisObs <- probAlive * probCapture
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture)
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen/probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) return(logProbData)
    return(exp(logProbData))
  }
)


## -----------------------------------------------------------------------------
C_dipper_model <- compileNimble(dipper_model)


## -----------------------------------------------------------------------------
C_dipper_model$calculate()
C_dipper_model$phi <- 0.7
C_dipper_model$p <- 0.5
C_dipper_model$calculate() # Ensure any lifted nodes are calculated
C_dipper_model$calculate('y[5,]')


## -----------------------------------------------------------------------------
dipper_MCMC <- buildMCMC(dipper_model)
C_dipper_MCMC <- compileNimble(dipper_MCMC, project = dipper_model)
samples <- runMCMC(C_dipper_MCMC)
summary(samples)


## -----------------------------------------------------------------------------
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma  # calculate once
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho) # vectorized alternative is given later
    return(result)
  })


## ---- include=FALSE-----------------------------------------------------------
# only needed for Rmd compilation; not needed for regular usage.
assign('expcov', expcov, .GlobalEnv)


## -----------------------------------------------------------------------------
code <- nimbleCode({
  mu[1:N] <- mu0 * ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # other parts of model omitted
})


## -----------------------------------------------------------------------------
code <- nimbleCode({
  # (hyper)parameter priors
  mu0 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  rho ~ dunif(0, 5)      # there might be a better non-informative prior for this

  # MVN normal (Gaussian process) prior
  mu[1:N] <- mu0 * ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  
  # likelihood for count data (e.g., disease mapping)
  for(i in 1:N) {
    expected[i] <- 1 # Data and/or predictive components would go here
    lambda[i] <- expected[i] * exp(x[i])
    y[i] ~ dpois(lambda[i])
  }
})

N <- 134
dists <- as.matrix(dist(runif(N)))
model <- nimbleModel(code, constants = list(N = N, dists = dists, ones = rep(1, N)), 
                     inits = list(rho = 1, sigma = 1, mu0 = 0))
deps <- model$getDependencies(c('rho','mu0','sigma'), self = FALSE)
deps
model$simulate(deps)  # may be a bit slow uncompiled given the nested looping
 # Note: there are no y values yet, so we are only looking at x
range(model$x)
model$calculate("x")
Cmodel <- compileNimble(model)
Cmodel$calculate("x")


## -----------------------------------------------------------------------------
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    result <- sigma*sigma * exp(-dists / rho)
    return(result)
  })

