#
#Joint analysis of three diseases
#
#1.- This script creates the data in a multivariate structure for INLA.
#2.- These data is tests using model with INLA.
#3.- The multivariate joint model described in the paper is fit with WinBUGS.

library("rgeos")
library("R2WinBUGS")
library("spdep")
library("maptools")
library("rgdal")
library("classInt")
library("RColorBrewer")
library(INLA)

#Load data
load("dismap_sim_data.RData")

#Rename data
spdf <- OyE.sim
rm(OyE.sim)

#Check names
names(spdf)
#[1] "obs.OCavity"   "esp.OCavity"   "obs.Esophagus" "esp.Esophagus"
#[5] "obs.Stomach"   "esp.Stomach"  

#Set shorter names
names(spdf) <- c("Obs.Cb", "Esp.Cb", "Obs.Eso", "Esp.Eso", 
  "Obs.Est", "Esp.Est")

#Create dataset for INLA (n x 3)
n <- nrow(spdf)
d <- list(OBS = matrix(NA, nrow = n*3, ncol = 3))
#Add observed
d$OBS[1:n, 1] <- spdf$Obs.Cb#Bucal cancer
d$OBS[n + 1:n, 2] <- spdf$Obs.Eso#Esophagous cancer
d$OBS[2*n + 1:n, 3] <- spdf$Obs.Est#Stomach cancer

#Expected cases
d$EXP <- c(spdf$Esp.Cb, spdf$Esp.Eso, spdf$Esp.Est)

#Area ID
d$AREAID <- rep(1:n, 3)

#Adjacency matrix
nb.spain <- poly2nb(spdf)
W <- nb2mat(nb.spain, style = "B")

#Define groups in data
d$r <- rep(1:3, each = n)
d$rf <- as.factor(d$r)

# Test some models with several spatial components with INLA 
# Note: THESE ARE NOT THE MODELS YOU ARE LOOKING FOR

#Formulas for models
form <- OBS ~ -1 + rf + f(AREAID, model = "besag", graph = W)
form.rep <- OBS ~ -1 + rf + f(AREAID, model = "besag", graph = W, replicate = r)
form.rep2 <- OBS ~ -1 + rf + f(AREAID, model = "besag", graph = W, replicate = r) + f(AREAID.com, delta, model = "besag", graph = W)

#Fit model
res1 <- inla(form, data = d, family = rep("poisson", 3), E = d$EXP)

#Fit model with replicated spatial effects
res2 <- inla(form.rep, data = d, family = rep("poisson", 3), E = d$EXP)


##
#Fit model with specific replicated spatial effects AND common spatial effects
##
#
#
#mlik for delta = 1 : -1067.407, -1067.925


#Define delta as the weighting for the common effect
#d$delta <- rep(1, 3*n)
delta.vec <- c(1.5, 2, 0.1) #Different loadings of the common effect
d$delta <- delta.vec[d$r]


d$AREAID.com <- d$AREAID
res3 <- inla(form.rep2, data = d, family = rep("poisson", 3), E = d$EXP,
  control.compute = list(dic = TRUE))

#Display summary results
spdf$SP.COM <- res3$summary.random$AREAID.com[, "mode"]
spdf$SP.BC <- res3$summary.random$AREAID[1:n, "mode"]
spdf$SP.EC <- res3$summary.random$AREAID[n + 1:n, "mode"]
spdf$SP.SC <- res3$summary.random$AREAID[2*n + 1:n, "mode"]

spplot(spdf, c("SP.BC", "SP.EC", "SP.SC", "SP.COM"))

#
#Fit the joint model using WinBUGS (this is the model described in the paper)
#
library(R2WinBUGS)
bugs.d <- list(obs = matrix(c(d$OBS[1:n, 1], d$OBS[1:n + n, 2], 
  d$OBS[1:n + 2*n, 3]), byrow = FALSE, ncol = 3))
bugs.d$exp <-  matrix(d$EXP, byrow = FALSE, ncol = 3)
bugs.d$D <- 3
bugs.d$N <- nrow(bugs.d$obs)
#ADD ADJACENCY MATRIX
bugs.d <- c(bugs.d, nb2WB(nb.spain))


#Initial values
bugs.inits <- list(alpha = rep(0, bugs.d$D), delta = rep(1, bugs.d$D),
  v = rep(0, bugs.d$N), spec = t(matrix(0, ncol = bugs.d$D, nrow = bugs.d$N)),
  tau.v = 1, tau.spec = 1)

BugsDir <- paste0(Sys.getenv("HOME"), 
  "/.wine/dosdevices/c:/Program Files/WinBUGS14")
model.file <- paste0(getwd(), "/joint_model.bug")
MCMCres<- bugs(data = bugs.d, inits = list(bugs.inits),
   working.directory = paste0(getwd(), "/WB"),
   parameters.to.save = c("alpha", "v", "spec", "delta", "tau.v", "tau.spec"),
   n.chains = 1, n.iter = 210000, n.burnin = 10000, n.thin = 20,
   model.file = model.file,
   bugs.directory = BugsDir,
   #debug = TRUE,
   WINEPATH="/usr/local/bin/winepath")


#Add spatial effects
spdf$SP.COMWB <- MCMCres$mean$v
spdf$SP.BCWB <- MCMCres$mean$spec[1, ]
spdf$SP.ECWB <- MCMCres$mean$spec[2, ]
spdf$SP.SCWB <- MCMCres$mean$spec[3, ]

#Display estimates of delta
par(mfrow = c(1,3))
for(i in 1:3) {
  plot(density(MCMCres$sims.list$delta[, i]), main = paste0("delta_",i))
}

spplot(spdf, c("SP.BC", "SP.EC", "SP.SC", "SP.COM", 
  "SP.BCWB", "SP.ECWB", "SP.SCWB", "SP.COMWB"))

#Save all results. These are used by other R files.
save(file = "joint.RData", list = ls())

