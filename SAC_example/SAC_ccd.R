#Example of use of CCD for integrating the hyperparamters

#Load packags
library(INLA)
library(INLABMA)
library(rsm)
library(parallel)

#Load fit.inla and other stuffprevious data
load("INLA-MH-SAC_data.RData")

#Obtain posterior modes from MCMC samples
#rho <- 0.34
#sd.rho <- 0.22
#lambda <- 0.32
#sd.lambda <- 0.28

#Values taken from sacsarlm output
#NOTE: These are actually ML estimates and NOT the modes
rho <- 0.35326
sd.rho <- 0.19669
lambda <- 0.13199
sd.lambda <- 0.29905

#Create CCD for spatial autocorrelation

w1 <-0.75 
w2 <- 0.75
#x1: lambda (error term)
#x2: rho (response term)
ccd.str <- ccd(2, n0 = 1, alpha = "rotatable", inscribed = TRUE,
  coding = list (
    x1 ~ (rho - rho)/(w1 * sd.rho), x2 ~ (lambda - lambda)/(w1 * sd.lambda))
)
ccd.str <- decode.data(ccd.str)

#Plot CCD
plot(ccd.str[, 3:4], pch = 19, xlim = c(-1, 1), ylim = c(-1,1))


#Fit models in parallel
options(mc.cores = 4)

#The output is similar to that one from INLAMH but with fewer models
res.ccd <- mclapply(1:nrow(ccd.str), function(X) {
  d <- fit.inla(d, unlist(as.vector(ccd.str[X, 3:4])))

  #Improve estimates of hyperparameters
  d$model <- inla.hyperpar(d$model)
  d$mlik <- d$model$mlik[1,1]
  d
})

#Marginal log-likelihood
mliks.ccd <- unlist(as.vector(lapply(res.ccd, function(X){X$mlik})))

#Posterior probabilities at CDD points
#The prior is a Uniform and it is not added
probs.ccd <- exp(mliks.ccd - max(mliks.ccd)) #Flat prior on rho and lambda
probs.ccd <- unlist(as.vector(probs.ccd/sum(probs.ccd)))

#Post means
sum(ccd.str$rho  * probs.ccd)
sum(ccd.str$lambda * probs.ccd)


#BMA to estimate the posterior marginals of the other parameters in the model
library(INLABMA)

#Fitted models
models.ccd <- lapply(res.ccd, function(X){X$model})
#Weights for BMA
ws.ccd <- probs.ccd

listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff.ccd <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models.ccd, ws.ccd, X)
    })

#Summary estimates
do.call(rbind, lapply(margeff.ccd[[1]], inla.zmarginal))

#Compute post. marginal of variance
marg.var.ccd <- inla.tmarginal(function(x) {1/x}, margeff.ccd[[2]][[1]])
inla.zmarginal(marg.var.ccd)

#Plot marginals
par(mfrow = c(2, 3))
lapply(margeff.ccd[[1]], function(X){plot(X, type = "l")})

#Variance
plot(marg.var.ccd, type = "l")
plot(density(ccd.str$rho, weights = probs.ccd, bw = .25) )
plot(density(ccd.str$lambda, weights = probs.ccd, bw = .25) )



#Compute impacts
rholambda.ccd <- ccd.str[, c("rho", "lambda")]

#Impacts
library(parallel)
options(mc.cores = 3)

#Total impacts: beta_r * (1/(1-rho))
aux <- sapply(rholambda.ccd[,1], function(X){1/(1-X)})

n.models <- length(models.ccd)
totimp.ccd <- mclapply(1:n.models, function(X) {
    aux[X] * models.ccd[[X]]$summary.fixed[-1, 1:2]
})
totimp.ccd <- Reduce("+", totimp.ccd) /n.models

#Direct impacts: trace ([I - rho * W]^{-1}) * beta_r
aux2 <- unlist(mclapply(rholambda.ccd[,1], function(X){
    mean(diag(solve(Diagonal(nrow(W)) - X * W)))
}))
dirimp.ccd <- mclapply(1:n.models, function(X) {
    aux2[X] * models.ccd[[X]]$summary.fixed[-1, 1:2]
})
dirimp.ccd <- Reduce("+", dirimp.ccd) / n.models

#Indirect impacts
indimp.ccd <- totimp.ccd - dirimp.ccd
indimp.ccd[, 2] <- sqrt(totimp.ccd[, 2]^2 - dirimp.ccd[ ,2]^2)

dirimp.ccd
indimp.ccd
totimp.ccd


save(file = "SAC_ccd.RData", list = ls())
