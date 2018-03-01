# Manski's model (SAC model) for spatial econometrics

The R code in this folder fits the SAC model to the Columbus dataset 
using MCMC, INLA within MCMC, INLA with CCD and INLA at maximum likelihood
estimates of the spatial autocorrelation parameters.

In addition to the R-INLA package, packages `INLABMA` and `SEMCMC` will be
required. You can install the latest versions as follows:

```
install.packages("INLA", repos=c(getOption("repos"), 
  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("INLABMA", repos="http://R-Forge.R-project.org")
devtools::install_github("becarioprecario/SEMCMC")
``` 

The main files included are:

* `INLA-MH-SAC.R`, implements INLA within MCMC and MCMC.
* `SAC_ccd.R`, implements INLA with CCD.
* `SAC_ml.R`, implements INLA at ML.
* `compute_impacts.R`, computes impacts for INLA within MCMC and MCMC.
* `plots.R`, creates some plots.

Files should be run in the order listed above. Also, `INLA-MH-SAC.R`
has been set with a lower number of iterations to run INLA within MCMC
than the ones used to produce the results in the paper. See


References:

Gómez-Rubio and Palmí-Perales (2017). Spatial Models with the Integrated Nested Laplace Approximation within Markov Chain Monte Carlo. [ArXiv paper arXiv:1702.03891](https://arxiv.org/abs/1702.03891). 
