#Implementation of INLA within M-H for Bayesian Inference
#  SAC MODEL (see INLABMA::sem.inla)

#Load libraries
library(INLA)

#This may be required to run the code on a cluster
#  INLA:::inla.dynload.workaround()

library(spdep)
library(INLABMA)#Includes INLAMH

#Load data and compute adjacency matrix
data(columbus)
lw <- nb2listw(col.gal.nb, style="W")

#FIT model using ML
sacml <- sacsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
  quiet=FALSE)


#Adjacency matrix
W <- as(as_dgRMatrix_listw(nb2listw(col.gal.nb)), "CsparseMatrix")
#Index for spatial random effects
columbus$idx<-1:nrow(columbus)
     
#Formula
form<- CRIME ~ INC + HOVAL
     
zero.variance = list(prec=list(initial = 25, fixed=TRUE))


#Fit model for a given rholambda = c(rho, lambda)
fit.inla <- function(data, rholambda) {
  res <- sac.inla(form, d = data, W.rho = W, W.lambda = W, 
    rho = rholambda[1], lambda = rholambda[2],
    family = "gaussian", impacts = FALSE,
    control.family = list(hyper = zero.variance),
    verbose = FALSE)

   return(list(mlik = res$mlik[1,1], model = res))
}


#Proposal distribution x -> y
#density
dq.rholambda <- function(x, y, sigma = .25, log =TRUE) {
	res <- dnorm(y, mean = x, sd = sigma, log = log)
  if(log) {
    return(sum(res))
  } else {

    return(prod(res))
  }
}
#random
rq.rholambda <- function(x, sigma = .25) {
	rnorm(length(x), mean = x, sd = sigma)
}

#Prior for beta: Uniform [-1.5, 1] (based on the eigenvalues of W)
prior.rholambda <- function(x, log = TRUE) {
  res <- dunif(x, -1.5, 1, log = log)
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

#Data as 'd'
d <- columbus

#Save data and functions so far. These are used in other R files.
save(file = "INLA-MH-SAC_data.RData", list = ls())

#Run INLA within MCMC: Settings used in the paper
#
#WARNING: This can take a long time and you may want to reduce the
#  number of iterations for testing purposes.
#inlamh.res <- INLAMH(d, fit.inla, c(0, 0), rq.rholambda, dq.rholambda, 
#  prior.rholambda, 
#  n.sim = 10000, n.burnin = 500, n.thin = 10, verbose = TRUE)

#Use this setting for testing:
# 1.- Starts at posterior modes (so convergence is faster)
# 2.- Uses less samples
# 3.- Some results may differ from the paper
inlamh.res <- INLAMH(d, fit.inla, c(0.4, 0.2), rq.rholambda, dq.rholambda, 
  prior.rholambda, 
  n.sim = 200, n.burnin = 20, n.thin = 1, verbose = TRUE)


#Save results just in case something fails later
#WARNING: This will create a HUGE file.
save(file = "INLA-MH-SAC.RData", list = ls())

#Display acceptance rate
print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 

#Attach results to access models and sampled values of c(rho, lambda)
attach(inlamh.res) 

#Obtain sampled values of c(rho, lambda)
rholambda.sim <- data.frame(do.call(rbind, b.sim))
names(rholambda.sim) <- c("rho", "lambda")

#Fit model using MCMC using packages rjags and SEMCMC
library(SEMCMC)
if(!file.exists("jags.RData")) {
  sac.mcmc <- SEMCMC(form, data = d, W = as.matrix(W), model = "sac",
    n.burnin = 500, n.iter = 10000, n.thin = 10)
  save(file = "jags.RData", list = "sac.mcmc")
} else {
  load("jags.RData")
}


#Plot posterior marginals of rho and lambda
pdf(file = "rho-lambda.pdf", , width = 10, height = 5)
  par(mfrow = c(1, 2))
  #rho
  plot(density(rholambda.sim$rho, bw = 0.1), 
    main = expression(rho), xlab = expression(rho))
  lines(density(unlist((sac.mcmc$results)[, "rho"]), bw = 0.1),
    col = "black", lty = 2)
  abline( v = coef(sacml)["rho"], lty = 4, lwd = 2)
  legend("topleft", legend = c("INLA w/ MCMC", "MCMC", "max. lik."),
    lty = c(1, 2, 4), lwd = c(1, 1, 2), 
    col = c("black", "black", "black"), bty = "n")

  #lambda
  plot(density(rholambda.sim$lambda, bw = 0.1), 
    main = expression(lambda), xlab = expression(lambda))
  lines(density(unlist((sac.mcmc$results)[, "lambda"]), bw = 0.1),
    col = "black", lty = 2)
  abline( v = coef(sacml)["lambda"], lty = 4, lwd = 2)
  legend("topleft", legend = c("INLA w/ MCMC", "MCMC", "max. lik."),
    lty = c(1, 2, 4), lwd = c(1, 1, 2),
    col = c("black", "black", "black"), bty = "n")
dev.off()

#Estimate joint posterior distribution using 2D kernel smoothing
library(MASS)
z.inla <- kde2d(rholambda.sim[,1], rholambda.sim[, 2])
z.mcmc <- kde2d(unlist((sac.mcmc$results)[, "rho"]),
  unlist((sac.mcmc$results)[, "lambda"]))

pdf(file = "contour.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))
  contour(z.inla, xlab = expression(rho), 
    ylab = expression(lambda), xlim = c(-.4, .8),
    ylim = c(-.4, .8), zlim = c(0, 5), main = "INLA within MCMC")
  points(coef(sacml)["rho"], coef(sacml)["lambda"], pch = 19)
  contour(z.mcmc, xlab = expression(rho), ylab = expression(lambda),
    lty = 2, xlim = c(-.4, .8),
    ylim = c(-.4, .8), zlim = c(0, 5), main = "MCMC")
  points(coef(sacml)["rho"], coef(sacml)["lambda"], pch = 19)
dev.off()


#Compute the posterior marginals of the other parameters in the model
#with Bayesian model averaging (BMA)
library(INLABMA)

#Obtain models and weights for BMA
models <- lapply(model.sim, function(X){X$model})
ws <- rep(1/length(models), length(models))

#BMA for posterior marginals of fixed effects and hyperparameters
listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models, ws, X)
    })

#Summary estimates
do.call(rbind, lapply(margeff[[1]], inla.zmarginal))

#Display posterior marginals of other parameters
pdf(file = "beta-tau.pdf")
  par(mfrow = c(2, 2))

  var.names <- c("Intercept", attr(terms(form), "term.labels"))
  #Coefficients
  for(i in 1:3) {
    plot(margeff[[1]][[i]], type = "l", main = var.names[i],
     xlab = "", ylab = "Density")
    lines(density(unlist((sac.mcmc$results)[, i])), col = "black", lty = 2)
    abline(v = coef(sacml)[i+2], lty = 3)
    legend("topleft", legend = c("INLA w/ MCMC", "MCMC", "max. lik."),
      lty = c(1, 2, 3), col = c("black", "black", "black"), bty = "n",
      cex = .8)
  }
  #Variance
  var.marg <- inla.tmarginal(function(x) {1/x}, margeff[[2]][[1]]) 
  plot(var.marg, type = "l", main = "Variance", xlab = "",
  ylab = "Density")
  lines(density(1/unlist((sac.mcmc$results)[, "tau"])), col = "black", lty = 2)
  abline(v = sacml$s2, lty = 3)
  legend("topright", legend = c("INLA w/ MCMC", "MCMC", "max. lik."),
    lty = c(1, 2, 3), col = c("black", "black", "black"), bty = "n",
    cex = .8)
dev.off()

#Save all results
#  THis is a HUGE file
save(file = "INLA-MH-SAC.RData", list = ls())

#Save summary results
rm(inlamh.res)
rm(models)
save(file = "INLA-MH-SAC_summary.RData", list = ls())
