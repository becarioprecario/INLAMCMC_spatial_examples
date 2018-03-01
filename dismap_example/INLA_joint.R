#
#Implementation of Bayesian shared components model with R-INLA and MCMC
#
#The variable estimated using MCMC is the 'delta' parameter
#
#

#Load libraries
library(sp)
library(spdep)
library(INLABMA)
library(INLA)
#This may be required to run the R script on a cluster
#  INLA:::inla.dynload.workaround()

library(R2WinBUGS)

#Load data from joint.R
#This contains the data structures required to fit the model with INLA
load("joint.RData")


#Fit linear model with R-INLA for a fixed 'delta'
#IMPORTANT: Only an intrinsic CAR is implemented in the shared component
#data: data used for model fitting: O, E, area, idx1 (1...0), idx2 (0...1)
#m.delta: value of delta
fit.inla <- function(data, m.delta) {
  data$delta <- m.delta[data$r]
  form.joint <- OBS ~ -1 + rf + 
    f(AREAID, model = "besag", graph = W, replicate = r,
      hyper = list(prec = list(param = c(0.01, 0.01)))) + 
    f(AREAID.com, delta, model = "besag", graph = W,
      hyper = list(prec = list(param = c(0.01, 0.01)))) 

  res <- inla(form.joint, data = data, family = rep("poisson", 3), E = data$EXP)

  return(list(mlik = res$mlik[1,1], model = res))

}


#Proposal distribution x -> y
#density function
dq.delta <- function(x, y, sigma = .5, log = TRUE) {
  res <- dlnorm(y, meanlog = log(x), sdlog = sigma, log = log)	

  if(log) {
    return(sum(res)) 
  } else {
    return(prod(res))
  }
}
#Generate random values
rq.delta <- function(x, sigma = .5) {
  rlnorm(length(x), meanlog = log(x), sdlog = sigma)
}

#Prior for delta (log-normal)
prior.delta <- function(x, mulog = 0, sigmalog = sqrt(1/.1), log = TRUE) {
  res <- dlnorm(x, meanlog = mulog, sdlog = sigmalog, log = log)

 if(log) {
    return(sum(res)) 
  } else {
    return(prod(res))
  }
}

#Fit model using INLA within MCMC: setting used in paper
# This TAKES A LONG
#inlamh.res <- INLAMH(d, fit.inla, c(1, 1, 1), rq.delta, dq.delta, prior.delta,
#  n.sim = 10000, n.burnin = 500, n.thin = 10, verbose = TRUE)

#Reduced number of simulations for testing
# Values start at posterio modes (approx.)
inlamh.res <- INLAMH(d, fit.inla, c(0.2522, 0.2753, 0.0558), rq.delta, 
  dq.delta, prior.delta,
  n.sim = 100, n.burnin = 20, n.thin = 1, verbose = TRUE)


#Save results just in case something goes wrong later
#Note: This will create a HUGE file
save(file = "INLA_joint.RData", list = ls())


#Create summary statistics and compute posterior marginals
attach(inlamh.res)

delta.sim2 <- do.call(rbind, b.sim)

#Summary of MCMC output and post. marginals of delta_i
par(mfcol = c(2,3))
for(i in 1:3) {
  delta.lab <-  bquote(delta[.(i)])
  plot(density(delta.sim2[, i]), main = delta.lab,
    xlab = delta.lab)
  lines(density(MCMCres$sims.list$delta[,i]), col = "red")
  abline(v = 1, lty = 3)

  plot(log(delta.sim2[, i]), type = "l", main = delta.lab, xlab = delta.lab)
}

#Display posterior marginals of weights delta_i
pdf(file = "delta.pdf", width = 10, height = 5)
bws1 <- c(0.25, 0.25, 0.25)
bws2 <- c(0.25, 0.25, 0.25)
par(mfcol = c(1,3))
for(i in 1:3) {
  delta.lab <-  bquote(delta[.(i)])
  plot(density(delta.sim2[, i], bw = bws1[i]), 
    main = delta.lab, xlab = delta.lab)
  lines(density(MCMCres$sims.list$delta[,i], bw = bws2[i]), lty = 3)#col = "red")
}
dev.off()


#BMA models to compute post. marginals of all the other model parameters
library(INLABMA)


#Fitted models with INLA
models <- lapply(model.sim, function(X){X$model})
#Weights for BMA
ws <- rep(1/length(models), length(models))

#BMA of fixed effects and hyperparameters to obtain their post. marginals
listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models, ws, X)
    })

#Summary estimates form post. marginals
do.call(rbind, lapply(margeff[[1]], inla.zmarginal))
do.call(rbind, lapply(margeff[[2]], inla.zmarginal))


#Marginals of the other parameters in the model
pdf(file = "params.pdf")
#Display fixed effects
par(mfrow = c(2,3))
for(i in 1:3) {
  plot(margeff[[1]][[i]], main = bquote(alpha[.(i)]), xlab = "", type = "l",
    ylab = "Density")
  lines(density(MCMCres$sims.list$alpha[, i]), lty = 3)#col = "red")
}

#Display averaged precision
#par(mfrow = c(1,2))
#Specific effect
plot(margeff[[2]][[1]], type = "l", xlab = "", main = "Prec. spec. effs.",
    ylab = "Density")
lines(density(MCMCres$sims.list$tau.spec), lty = 3)#col = "red")

#Common effect
plot(margeff[[2]][[2]], type = "l", xlab = "", main = "Prec. common eff.", 
  xlim = c(0, 100), ylab = "Density")
lines(density(MCMCres$sims.list$tau.v), lty = 3)#col = "red")
dev.off()

#Save all results
#Note: This is a HUGE file.
save(file = "INLA_joint.RData", list = ls())

#Keep only summary stuff
rm(inlamh.res)
rm(models)
save(file = "INLA_joint_summary.RData", list = ls())


