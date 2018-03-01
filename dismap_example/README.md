# Multivariate disease mapping model with shared components

The R code in this folder fits the multivariate disease mapping model 
with a shared spatial random effect and disease-specific spatial random effects
to a *simulated* dataset in Spain at the province level.
Estimations methods are
MCMC, INLA within MCMC, INLA with CCD and INLA at the posterior modes
of the disease-specific weights $\delta^{(d)}$.

In addition to the R-INLA package, package `INLABMA` will be
required. You can install the latest versions as follows:

```
install.packages("INLA", repos=c(getOption("repos"), 
  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("INLABMA", repos="http://R-Forge.R-project.org")
``` 

Files included are:

* `dismap_sim_data.RData`, simulated data. Observed and expected cases
of oral cavity, esophagus and stomach cancer in Spain in the period 1996-2014 at the province level.
* `joint.R`, load data and create data structure for INLA and fit model with WinBUGS.
* `INLA_joint.R`, implementes INLA within MCMC.
* `joint_ccd.R`, implements INLA with CCD.
* `joint_mode.R`, implements INLA at MODE.
* `plots.R`, produces some plots of the posterior marginals.
* `maps.R`, create some maps with summary information.



References:

Gómez-Rubio and Palmí-Perales (2017). Spatial Models with the Integrated Nested Laplace Approximation within Markov Chain Monte Carlo. [ArXiv paper arXiv:1702.03891](https://arxiv.org/abs/1702.03891).

