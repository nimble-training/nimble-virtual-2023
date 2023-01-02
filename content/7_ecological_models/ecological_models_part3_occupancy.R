## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=FALSE)
library(nimble)


## ---- eval=FALSE--------------------------------------------------------------
## y[i, 1:T] ~ dOcc_v(probOcc = psi[i],
##                    probDetect = p[i, 1:T], len = T)


## ---- eval=FALSE--------------------------------------------------------------
## y[1:num_sites, 1:num_years] ~ dDynOcc_svv(probPersist = psi,
##                                           probColonize = gamma[1:num_years],
##                                           p = p[1:num_years])

