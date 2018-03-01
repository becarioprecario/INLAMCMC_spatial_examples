#
#Make some maps with the output
#

library(INLA)
library(R2WinBUGS)
library(spdep)

#Spatial common effects
load("INLA_joint.RData")

#Load data
load("joint.RData")

#Get models
models <- lapply(inlamh.res$model.sim, function(X){X$model})


#SMR
spdf$SMR.Cb <- spdf$Obs.Cb/spdf$Esp.Cb
spdf$SMR.Eso <- spdf$Obs.Eso/spdf$Esp.Eso
spdf$SMR.Est <- spdf$Obs.Est/spdf$Esp.Est


pdf(file = "SMR.pdf", width = 10, height = 5)
n.cuts <- 9
spplot(spdf, c("SMR.Cb", "SMR.Eso", "SMR.Est"), 
  names.attr = c("Lip, Oral Cavity & Pharynx Cancer", "Esophagus Cancer", "Stomach Cancer"),
  cuts = n.cuts, col.regions = rev(grey.colors(n.cuts + 1)))
dev.off()

#
#Summarise posterior median of shared spatial random effect
#

#INLA within MCMC
aux <- sapply(1:length(models), function (X) {
  models[[X]]$summary.random$AREAID.com[, "0.5quant"]
})
aux <- apply(aux, 1, median)

spdf@data$COMMON.INLAMCMC <- aux
#MCMC
spdf@data$COMMON.WB <- MCMCres$median$v

#INLA with CCD
load("joint_ccd.RData")
spdf@data$COMMON.INLACCD <- sapply(1:47, function(i) {
  margs <- lapply(models.ccd, function(X){X$marginals.random[[2]][[i]]})
  inla.zmarginal(INLABMA:::fitmargBMA(margs, ws), TRUE)$quant0.5
})

#INLA at mode
load("joint_mode.RData")
spdf@data$COMMON.INLAMODE <- model.mode$model$summary.random[[2]][, "0.5quant"]

#Display posterior median of shared spatial random effect
pdf(file = "common_all.pdf", width = 10, height = 10)
  n.cuts <- 19
  spplot(spdf, c("COMMON.INLAMCMC", "COMMON.INLACCD", "COMMON.INLAMODE", 
    "COMMON.WB"), 
    names.attr = c("INLA w/ MCMC", "INLA w/ CCD", "INLA at MODE", "MCMC"),
    cuts = n.cuts, col.regions = rev(grey.colors(n.cuts + 1)))
dev.off()

