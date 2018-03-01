#Implementation of INLA within M-H for Bayesian Inference
#  SAC MODEL (see INLABMA::sem.inla)

#Load libraries
library(INLA)

#To be used when running on cluster
#  INLA:::inla.dynload.workaround()


library(spdep)
library(INLABMA)#Includes INLAMH

load(file = "INLA-MH-SAC.RData")
load("jags.RData")
#load("impacts.RData") #Load previous impacts; NOT TO compute MCMC impacts

print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 

#Reduce size for testing
#inlamh.res$b.sim <- inlamh.res$b.sim[1:100]
#inlamh.res$model.sim <- inlamh.res$model.sim[1:100]

attach(inlamh.res) 

rholambda.sim <- data.frame(do.call(rbind, b.sim))#[1:9998,]
names(rholambda.sim) <- c("rho", "lambda")

#Impacts
library(parallel)
options(mc.cores = 4)

#Remove weird modes fitted
mliks <- unlist(lapply(model.sim, function(X){X$mlik}))
idx <- which(mliks < 0)

#Total impacts: beta_r * (1/(1-rho))
aux <- sapply(rholambda.sim[,1], function(X){1/(1-X)})

n.models <- length(model.sim[idx])
totimp <- mclapply(idx, function(X) {
    aux[X] * cbind(model.sim[[X]]$model$summary.fixed[-1, 1],
      aux[X] * model.sim[[X]]$model$summary.fixed[-1, 2]^2)
})
totimp <- Reduce("+", totimp) /n.models

#Direct impacts: trace ([I - rho * W]^{-1}) * beta_r
aux2 <- unlist(mclapply(rholambda.sim[,1], function(X){
    mean(diag(solve(Diagonal(nrow(W)) - X * W)))
}))
dirimp <- mclapply(idx, function(X) {
    aux2[X] * cbind(model.sim[[X]]$model$summary.fixed[-1, 1],
      aux2[X] * model.sim[[X]]$model$summary.fixed[-1, 2]^2)
})
dirimp <- Reduce("+", dirimp) / n.models

#Indirect impacts
indimp <- totimp - dirimp
indimp[, 2] <- totimp[, 2] - dirimp[ ,2]

#Compute standard deviations
dirimp[, 2] <- sqrt(dirimp[, 2])
indimp[, 2] <- sqrt(indimp[, 2])
totimp[, 2] <- sqrt(totimp[, 2])

dirimp
indimp
totimp

#Max. Lik.
impacts.maxlik <- spdep::impacts(sacml, listw = lw)
impacts.maxlik

#Compute impacts for MCMC model
#  This may take a while...
library(SEMCMC)
impacts.mcmc <- impacts(sac.mcmc, W)
impacts.mcmc


save(file = "impacts.RData",
  list = c("dirimp", "indimp", "totimp", "impacts.maxlik", "impacts.mcmc"))
