#Plots of posterior marginals and bivariate joint posterior distribution
#  of rho and lambda.

#Load libraries
library(INLA)

#To be used when running on cluster
#  INLA:::inla.dynload.workaround()


library(spdep)
library(INLABMA)#Includes INLAMH

#Load previous results
load("SAC_ccd.RData")
load("jags.RData")
load("SAC_ml.RData")
load("INLA-MH-SAC_summary.RData")

#Samples from MCMC
rho.mcmc <- unlist(sac.mcmc$results[,"rho"])
lambda.mcmc <- unlist(sac.mcmc$results[,"lambda"])

#Plot posterior marginals of rho and lambda
pdf(file = "rho-lambda_all.pdf", , width = 10, height = 5)
  par(mfrow = c(1, 2))
  #rho
  plot(density(rholambda.sim$rho, bw = 0.1), 
    main = expression(rho), xlab = expression(rho))
  lines(density(ccd.str$rho, weights = probs.ccd, bw = .25), lty = 2 )
  lines(density(rho.mcmc, bw = 0.1),
     col = "black", lty = 3)
  abline( v = coef(sacml)["rho"], lty = 4, lwd = 2)
  legend("topleft", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", "max. lik."),
    lty = c(1, 2, 3, 4), lwd = c(1, 1, 1, 2), 
    col = "black",  bty = "n")

  #lambda
  plot(density(rholambda.sim$lambda, bw = 0.1), 
    main = expression(lambda), xlab = expression(lambda))
  lines(density(ccd.str$lambda, weights = probs.ccd, bw = .25), lty = 2 )
  lines(density(lambda.mcmc, bw = 0.1),
    col = "black", lty = 3)
  abline( v = coef(sacml)["lambda"], lty = 4, lwd = 2)
  legend("topleft", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", "max. lik."),
    lty = c(1, 2, 3, 4), lwd = c(1, 1, 1, 2),
    col = "black", bty = "n")
dev.off()

#Display joint posterior distribution of rho and lambda
# using MCMC and INLA within MCMC
library(MASS)
z.inla <- kde2d(rholambda.sim[,1], rholambda.sim[, 2])
z.mcmc <- kde2d(rho.mcmc, lambda.mcmc)
  

pdf(file = "contour_all.pdf", width = 10, height = 5)
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


#Display posterior marginals of other parameters
pdf(file = "beta-tau_all.pdf")
  par(mfrow = c(2, 2))

  var.names <- c("Intercept", attr(terms(form), "term.labels"))
  #Coefficients
  for(i in 1:3) {
    plot(margeff[[1]][[i]], type = "l", main = var.names[i],
     xlab = "", ylab = "Density")
    lines(margeff.ccd[[1]][[i]], lty = 2)
    lines(density(unlist(sac.mcmc$results[, i])), col = "black", lty = 3)
    lines(model.ml$model$marginals.fixed[[i]], lty = 4)
    abline(v = coef(sacml)[i+2], lty = 4)
    legend("topleft", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", "INLA at ML"),
      lty = c(1, 2, 3, 4), col = rep("black", 4), bty = "n",
      cex = .8)
  }
  #Variance
  var.marg <- inla.tmarginal(function(x) {1/x}, margeff[[2]][[1]]) 
  plot(var.marg, type = "l", main = "Variance", xlab = "",
  ylab = "Density")
  lines(marg.var.ccd, lty = 2)
  lines(density(1/unlist(sac.mcmc$results[, "tau"])), col = "black", lty = 3)
  lines(inla.tmarginal(function(x) {1/x},
    model.ml$model$marginals.hyperpar[[1]]), col = "black", lty = 4)
  abline(v = sacml$s2, lty = 4)
  legend("topright", legend = c("INLA w/ MCMC", "INLA w/ CCD", "MCMC", "INLA at ML"),
    lty = c(1, 2, 3, 4), col = c("black", "black", "black"), bty = "n",
    cex = .8)
dev.off()
