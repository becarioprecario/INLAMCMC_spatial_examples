#
#Create plots with estimates from all four estimation methods described in the paper
#
#The variable estimated using MCMC is the 'delta' parameter
#

#Load libraries
library(sp)
library(spdep)
library(INLABMA)
library(INLA)
#To be used when running on cluster
#  INLA:::inla.dynload.workaround()

library(R2WinBUGS)

#Load data from joint.R
load("joint.RData")
load("joint_ccd.RData")
load("INLA_joint_summary.RData")#"INLA_joint.RData" minus 'inlamh.res' and 'models'
load("joint_mode.RData")

#Display post. marginals of delta with MCMC and INLA within MCMC
par(mfcol = c(2,3))
for(i in 1:3) {
  plot(density(delta.sim2[, i]), main = bquote(delta^(.(i))), xlab = "x")
  lines(density(MCMCres$sims.list$delta[,i]), col = "red")
  abline(v = 1, lty = 3)
}

#Display estiamtes of post. marginal of delta_i with all 4 methods
pdf(file = "delta_all.pdf", width = 10, height = 5)
bws1 <- c(0.2, 0.2, 0.25)
bws2 <- bws1
par(mfcol = c(1,3))
for(i in 1:3) {
  plot(density(delta.sim2[, i], bw = bws1[i]), 
    main = bquote(delta^(.(i))), xlab = "")
  lines(density(ccd.str[, 2 + i], weights = probs.ccd, bw =  bws1[i]), lty = 2)
  lines(density(MCMCres$sims.list$delta[,i], bw = bws2[i]), lty = 3)#col = "red")
  abline(v = delta.mode[i], lty = 4)

  legend("topright", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", 
    "INLA at MODE"), bty = "n", lty = 1:4, cex = 1.2)
}
dev.off()

#Display post. marginals of the other parameters with all 4 methods
pdf(file = "params_all.pdf")
#Display fixed effects
par(mfrow = c(2,3))
for(i in 1:3) {
  plot(margeff[[1]][[i]], main = bquote(alpha^(.(i))), type = "l",
    xlab = "", ylab = "Density")
  lines(margeff.ccd[[1]][[i]], lty = 2)
  lines(density(MCMCres$sims.list$alpha[, i]), lty = 3)#col = "red")
  lines(model.mode$model$marginals.fixed[[i]], lty = 4)

  legend("topright", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", 
    "INLA at post. mode"), bty = "n", lty = 1:4, cex = 0.5)
}

#Display averaged precision
#par(mfrow = c(1,2))
#Precision of specific effect
plot(margeff[[2]][[1]], type = "l", main = expression(tau[s]),
  xlab = "", ylab = "Density")
lines(margeff.ccd[[2]][[1]], lty = 2)
lines(density(MCMCres$sims.list$tau.spec), lty = 3)#col = "red")
lines(model.mode$model$marginals.hyperpar[[1]], lty = 4)

legend("topright", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", 
  "INLA at post. mode"), bty = "n", lty = 1:4, cex = 0.85)

#Precision of shared spatial random effect effect
plot(margeff[[2]][[2]], type = "l", main = expression(tau[v]),
  xlim = c(0, 100), ylim = c(0, 0.10),
  xlab = "", ylab = "Density")
lines(margeff.ccd[[2]][[2]], lty = 2)
lines(density(MCMCres$sims.list$tau.v), lty = 3)#col = "red")
lines(model.mode$model$marginals.hyperpar[[2]], lty = 4)
legend("topright", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC",  
  "INLA at post. mode"), bty = "n", lty = 1:4, cex = 0.85)

dev.off()

