#Example of use of CCD for integrating the hyperparamters

#Load packags
library(INLA)
library(rsm)
library(parallel)

#Load previous data
load("INLA_joint_summary.RData")#INLA_joint.RData minus inlamh.res and models 

#Posterior means and standard deviations
# inla.zmarginal(do.call(cbind, density(delta.sim2[, 1])[c("x", "y")]))

#Get post. mode from kernel density estimate of post. marginal
#x: Samples from MCMC
get.mode <- function(x) {
  xy <- do.call(cbind, density(x, bw = 0.2)[c("x", "y")])
  xy[which.max(xy[, 2]), 1]
}
apply(delta.sim2, 2, get.mode)

delta.mode <- apply(delta.sim2, 2, get.mode)
delta.mean <- apply(delta.sim2, 2, mean)
delta.sd <- apply(delta.sim2, 2, sd) # 0.6635271 0.7262760 0.1733759
#delta.mode <- c(1, 1, 0.5) #c(0.15, 0.15, 0.01)
#delta.sd <- c(1, 1, 0.5)

#Set CCD in log-scale
log.scale <- TRUE

#Change to log-scale
if(log.scale ) {
  delta.sd <- apply(log(delta.sim2), 2, sd)
  #Delta method
  #delta.sd <- sqrt(delta.sd^2 / delta.mean^2)
  delta.mode <- log(delta.mode)
}
#Create CCD for spatial autocorrelation
#x1: lambda (error term)
#x2: rho (response term)
ws.ccd <- rep(0.75, 3) #Weights for ccd
ccd.str <- ccd(3, n0 = 1, alpha = "rotatable", inscribed = TRUE,
  coding = list (
#    x1 ~ (rho - mode.rho)/sd.rho, x2 ~ (lambda - mode.lambda)/sd.lambda)
    x1 ~ (delta.1 - delta.mode[1])/(ws.ccd[1] * delta.sd[1]),
    x2 ~ (delta.2 - delta.mode[2])/(ws.ccd[2] * delta.sd[2]),
    x3 ~ (delta.3 - delta.mode[3])/(ws.ccd[3] * delta.sd[3]))
)
ccd.str <- decode.data(ccd.str)

#Go back to original scale, if needed
if(log.scale) {
  ccd.str[, 3:5] <- exp(ccd.str[, 3:5])
} else {
  #If original scale is used, remove negative values of delta_i
  ccd.str <- ccd.str [with(ccd.str, delta.1 >0 & delta.2 > 0 & delta.3 > 0), ]
}
ccd.str


#Fit models in parallel at CCD points
options(mc.cores = 4)

res.ccd <- mclapply(1:nrow(ccd.str), function(X) {
  d <- fit.inla(d, unlist(as.vector(ccd.str[X, 3:5])))

  d$model <- inla.hyperpar(d$model)
  d$mlik <- d$model$mlik[1,1]
  d
})


#Marginal log-likelihood
mliks.ccd <- unlist(as.vector(lapply(res.ccd, function(X){X$mlik})))

#Compute posterior prob. at CCD points; these are also the weights for BMA
probs.ccd <- exp(mliks.ccd - max(mliks.ccd))
#Add prior
probs.ccd <- probs.ccd + exp(apply(ccd.str[, 3:5], 1,  prior.delta))
probs.ccd <- unlist(as.vector(probs.ccd/sum(probs.ccd)))

#Post means
sum(ccd.str$delta.1  * probs.ccd)
sum(ccd.str$delta.2 * probs.ccd)
sum(ccd.str$delta.3 * probs.ccd)


#BMA to compute the posterior marginals of the other parameters in the model
library(INLABMA)

#Models fitted with INLA
models.ccd <- lapply(res.ccd, function(X){X$model})
#Weights for BMA
ws <- probs.ccd

listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff.ccd <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models.ccd, ws, X)
    })

#xxr must be set to c(0, 100) manually for the second hyperpar in order to be
# estimated in a proper range.
#This must be fixed in INLABMA::fitmargBMA2
if(FALSE) {
  debug(INLABMA:::fitmargBMA2)
  kk <- INLABMA:::fitmargBMA2(models.ccd, ws, "marginals.hyperpar")
  margeff.ccd[2] <- list(kk)
}# if(FALSE)

#Summary estimates
do.call(rbind, lapply(margeff.ccd[[1]], inla.zmarginal))
#Variance
marg.vars.ccd <- lapply(margeff.ccd[[2]], function(X) {
  inla.tmarginal(function(x){ 1/x}, X)}
)
do.call(rbind, lapply(marg.vars.ccd, inla.zmarginal))

#Plot posterior marginals
par(mfrow = c(2, 3))
lapply(margeff.ccd[[1]], function(X){plot(X, type = "l")})
#Precisions
lapply(margeff.ccd[[2]], function(X){plot(X, type = "l")})

dev.new()
par(mfrow = c(1, 3))
plot(density(ccd.str$delta.1, weights = probs.ccd, bw = .35) )
plot(density(ccd.str$delta.2, weights = probs.ccd, bw = .35) )
plot(density(ccd.str$delta.3, weights = probs.ccd, bw = .25) )

save( file = "joint_ccd.RData", list = c("ccd.str", "models.ccd",
  "margeff.ccd", "mliks.ccd", "probs.ccd"))

#Summary statistics and correlations between delta_i's
means.ccd  <- apply(ccd.str[, 3:5], 2, function(X) {sum(X * probs.ccd)})
vars.ccd  <- apply(ccd.str[, 3:5], 2, function(X) {sum(X^2 * probs.ccd)}) -
   means.ccd^2

#Correlation
for(i in 1:2) {
  for(j in (i+1):3) {
    print(
    sum((ccd.str[, 2 + i] - means.ccd[i]) * (ccd.str[, 2 + j] - means.ccd[j]) *
      probs.ccd) / sqrt(vars.ccd[i] * vars.ccd[j])
    )
  }
}
