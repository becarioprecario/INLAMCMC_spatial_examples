#
#Fit model at mode of delta
#
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
load("INLA_joint_summary.RData")#INLA_joint.RData minus 'inlamh.res' and 'models'

#Fit model at mode
delta.mode <- c(0.2522, 0.2753, 0.0558)#joint_ccd.R; original dataset
model.mode <- fit.inla(d, delta.mode)
model.mode$model <- inla.hyperpar(model.mode$model)
model.mode$mlik <- model.mode$model$mlik[1,1]

save(file = "joint_mode.RData", list = c("delta.mode", "model.mode"))

