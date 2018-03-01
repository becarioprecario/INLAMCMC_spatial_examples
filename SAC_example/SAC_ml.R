#Fit SAC model using maximum likelihood (ML)
#  SAC MODEL (see INLABMA::sem.inla and spdep::sacsarlm)

#Load libraries
library(INLA)
library(INLABMA)

library(spdep)

#Load fit.inla, sacml and some other stuff
load("INLA-MH-SAC_data.RData")

#Load data and compute weight matrix
data(columbus)
lw <- nb2listw(col.gal.nb, style = "W")

#Fit model at ML estimates
ml.rho <- sacml$rho
ml.lambda <- sacml$lambda
model.ml <- fit.inla(d, c(ml.rho, ml.lambda))
#Improve estimation of hyperparameters
model.ml$model <- inla.hyperpar(model.ml$model)
model.ml$mlik <- model.ml$model$mlik[1,1]


#Compute Impacts
# rho and lambda are kept fixed.
library(parallel)
options(mc.cores = 4)

#Total impacts: beta_r * (1/(1-rho))
aux <- 1/(1 - ml.rho)

totimp.ml <- aux * model.ml$model$summary.fixed[-1, 1:2]

#Direct impacts: trace ([I - rho * W]^{-1}) * beta_r
aux2 <- mean(diag(solve(Diagonal(nrow(W)) - ml.rho * W)))
dirimp.ml <- aux2 * model.ml$model$summary.fixed[-1, 1:2]


#Indirect impacts
indimp.ml <- totimp.ml - dirimp.ml
indimp.ml[, 2] <- sqrt(totimp.ml[, 2]^2 - dirimp.ml[ ,2]^2)

#Show impacts
dirimp.ml
indimp.ml
totimp.ml


save(file = "SAC_ml.RData", list = c("model.ml", "dirimp.ml",
  "indimp.ml", "totimp.ml"))
