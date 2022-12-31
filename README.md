# nimble-virtual-2023

Materials for the virtual NIMBLE workshop, January 4-6, 2023. 

To prepare for the workshop:

 - Install NIMBLE (see below)
 - Install additional packages (see below)
 - Join our Slack workspace (highly recommended but not required; see information in email)
 - Download these materials (and check back before the workshop on Wednesday for updates)
 - See email for Zoom invitation.

All materials for the workshop will be in this GitHub repository. If you're familiar with Git/GitHub, you already know how to get all the materials on your computer. If you're not, simply click [here](https://github.com/nimble-training/nimble-virtual-2023/archive/main.zip).

There is some overview information [here (https://htmlpreview.github.io/?https://github.com/nimble-training/nimble-virtual-2023/blob/main/overview.html), including links to the content modules in order.

Time: 8 am - 1 pm California time (GMT-8); 11 am - 4 pm Eastern US time.

Location: Zoom (see email) and [Slack](https://2023nimbleworkshop.slack.com).

## Slack and Zoom instructions and tips/tricks

Please see [this document (https://docs.google.com/document/d/1hhm6Eco0KevM30aDGdpo0n-gYDcT7IPbw7_vTrb-BVQ/edit?usp=sharing) for discussion of how to use Zoom and Slack during the workshop.

## Tentative Schedule

All times are California time. If we need to modify the start times of
units, we'll try to announce in advance on the Slack #general channel.

Day 1 (Wednesday January 4):

1. (8 am - 8:45 am) Introduction to NIMBLE: Basic concepts and workflows
2. (9 am - 10 am) Working with NIMBLE models and converting from WinBUGS/JAGS
3. (10:30 am - 11:45 pm) Comparing and customizing MCMC methods in NIMBLE
4. (12 pm - 1 pm) Strategies for improving MCMC

Day 2 (Thursday January 5):

5. (8 am - 9 am) Writing your own functions and distributions 
6. (9:15 am - 10:30 am) Spatial modeling
7. (11 am - 1 pm) Ecological models (including Hidden Markov Models)

Day 3 (Friday January 6):

8. (8 am - 9 am) Introduction to nimbleFunctions and nimble programming
9. (9:15 am - 10 am) Introduction to automatic differentiation (AD) (with Laplace approximation example)
10. (10:30 am - 12 noon) nimbleFunction programming (with user-defined sampler example)
11. (12:15 pm - 1 pm) Deeper into programming with AD

## Help with NIMBLE

Our user manual is [here](https://r-nimble.org/html_manual/cha-welcome-nimble.html).

We have a 'cheatsheet' and a guide to converting from JAGS or WinBUGS to NIMBLE [here](https://r-nimble.org/documentation).

For the small number of you who are not already familiar with writing models in WinBUGS, JAGS, or NIMBLE, you may want to look through the first module (Introduction to NIMBLE) or Section 5.2 of our user manual in advance.

We're happy to answer questions about writing models as we proceed through the workshop, but if you have no experience with it, reviewing in advance will greatly lessen the odds you feel lost right at the beginning.

## Installing NIMBLE

NIMBLE is an R package on CRAN, so in general it will be straightforward to install as with any R package, but you do need a compiler and related tools on your system.  

In summary, here are the steps.

1. Install compiler tools on your system. [https://r-nimble.org/download](https://r-nimble.org/download) has more details on how to install *Rtools* on Windows and how to install the command line tools of *Xcode* on a Mac. Note that if you have packages requiring a compiler (e.g., *Rcpp*) on your computer, you should already have the compiler tools installed.

2. Install the *nimble* package from CRAN in the usual fashion for an R package. More details (including troubleshooting tips) can also be found in Section 4 of the [NIMBLE manual](https://r-nimble.org/html_manual/cha-installing-nimble.html).

3) To test that things are working please run the following code in R:

```
library(nimble)
code <- nimbleCode({
  y ~ dnorm(0,1)
})
model <- nimbleModel(code)
cModel <- compileNimble(model)
```


If that runs without error, you're all set. If not, please see the troubleshooting tips and email nimble.stats@gmail.com directly if you can't get things going.  

In general we encourage you to update to the most recent version of NIMBLE, 0.11.0.


#### (Required) Development version(s) of NIMBLE

We'll be using the new automatic differentiation (AD) features of NIMBLE.  These aren't yet in the version of NIMBLE released through CRAN, so you'll have to install it manually as follows in order to replicate the parts of the workshop that make use of AD, including HMC sampling and Laplace approximation.  

Please follow [these instructions](https://r-nimble.org/ad-beta) to install the beta version of NIMBLE with AD support, including the `nimbleHMC` package and a beta version of `nimbleEcology`.

## Installing additional packages

Some of the packages we will use (beyond those automatically installed with `nimble`) can be installed as follows:

```
install.packages(c("compareMCMCs", "mcmcplots", "CARBayesdata", "sp", "spdep", "classInt", "coda", "nimbleEcology", "nimbleSCR"))
```

