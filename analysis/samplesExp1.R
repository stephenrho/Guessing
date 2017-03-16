# run samples for experiment 1
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)
library(rjags)
library(R2jags)

# READ DATA AND MODELS ----------------------------------------------------
setwd('analysis/')

exp1 <- read.table(file = 'data/exp1/exp1_full')
exp1$Cond <- as.factor(exp1$Cond)

# STRUCTURE DATA EXP 1 -----------------------------------------------------
DATALIST_EXP1 <- list(
  y = exp1$respChange,
  CHANGE = exp1$testChange,
  C = rep(1, nrow(exp1)),
  B_level = as.numeric(exp1$Cond),
  B_n = length(unique(exp1$Cond)),
  N = exp1$SetSize,
  S = length(unique(exp1$ppt)),
  id = exp1$ppt,
  n = nrow(exp1),
  # for data model
  condLevel = as.numeric(exp1$Cond),
  N_c = length(unique(exp1$Cond)),
  ssLevel = as.numeric(as.factor(exp1$SetSize)),
  N_ss = length(unique(exp1$SetSize)),
  testLevel = as.numeric(as.factor(exp1$testChange))
)

# RUN CHAINS EXP 1 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD')

dataM.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(dataM.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD') 

inf.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/informedWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# variable k
inf_vk.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/informedWD_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf_vk.samplesE1$BUGSoutput$summary)$Rhat < 1.1)
# K[1] is unstable. Increasing iterations will not solve

# UNINFORMED
uninf.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/uninformedWD.txt", n.chains = 5, n.iter = 55000, n.burnin = 5000, n.thin = 5)
# additional iterations and thinning run to converge

all(as.data.frame(uninf.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# variable k
uninf_vk.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/uninformedWD_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf_vk.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'P_og_subj')

mixture.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/mixtureWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(mixture.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# logistic RULE
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Lsubj')

logistic.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/logisticWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(logistic.samplesE1$BUGSoutput$summary)$Rhat < 1.1)

# EXTRACT POSTERIOR --------------------------------------------------------
source('guessingFunctions.R')

dataM.mcmcE1 <- as.matrix(as.mcmc(dataM.samplesE1))
inf.mcmcE1 <- as.matrix(as.mcmc(inf.samplesE1))
uninf.mcmcE1 <- as.matrix(as.mcmc(uninf.samplesE1))
mixture.mcmcE1 <- as.matrix(as.mcmc(mixture.samplesE1))
logistic.mcmcE1 <- as.matrix(as.mcmc(logistic.samplesE1))

dataM.postE1 <- logistic(dataM.mcmcE1[,3:20])

# extract estimated median f,h pairs and HDIs
fh.estimates.postE1 <- expand.grid(N = c(2,5,8), P = c(0.3,0.5,0.7), 
                                   f.median=0, f.lower=0, f.upper=0,
                                   h.median=0, h.lower=0, h.upper=0)

for (p in 1:3){
  for (n in 1:3){
    fcol <- paste0('R[1,', p, ',', n, ']')
    hcol <- paste0('R[2,', p, ',', n, ']')
    fh.estimates.postE1[fh.estimates.postE1$N == c(2,5,8)[n] & fh.estimates.postE1$P == c(.3,.5,.7)[p],3:8] <- c(median(dataM.postE1[,fcol]), HDIofMCMC(dataM.postE1[,fcol]), median(dataM.postE1[,hcol]), HDIofMCMC(dataM.postE1[,hcol]))
  }
}

# INFORMED
# extract posterior on natural scale
inf.postE1 <- cbind(inf.mcmcE1[,'K'],
                    logistic(inf.mcmcE1[,'A']),
                    logistic(inf.mcmcE1[,'U[1]']),
                    logistic(inf.mcmcE1[,'U[2]']),
                    logistic(inf.mcmcE1[,'U[3]']))

# UNINFORMED
# extract posterior on natural scale
uninf.postE1 <- cbind(uninf.mcmcE1[,'K'],
                      logistic(uninf.mcmcE1[,'A']),
                      logistic(uninf.mcmcE1[,'U[1]']),
                      logistic(uninf.mcmcE1[,'U[2]']),
                      logistic(uninf.mcmcE1[,'U[3]']))

# MIXTURE
# extract posterior on natural scale
mixture.postE1 <- cbind(mixture.mcmcE1[,'K'],
                        logistic(mixture.mcmcE1[,'A']),
                        logistic(mixture.mcmcE1[,'U[1]']),
                        logistic(mixture.mcmcE1[,'U[2]']),
                        logistic(mixture.mcmcE1[,'U[3]']),
                        logistic(mixture.mcmcE1[,'POG']))

mixture_POG.postE1 <- logistic(mixture.mcmcE1[,paste0('P_og_subj[', 1:DATALIST_EXP1$S , ']')])
mixture_POG.postE1 <- cbind(apply(mixture_POG.postE1, MARGIN = 2, median), t(apply(mixture_POG.postE1, MARGIN = 2, HDIofMCMC)))

# logistic
logistic.postE1 <- cbind(logistic.mcmcE1[,'K'],
                        logistic(logistic.mcmcE1[,'A']),
                        logistic(logistic.mcmcE1[,'U[1]']),
                        logistic(logistic.mcmcE1[,'U[2]']),
                        logistic(logistic.mcmcE1[,'U[3]']),
                        exp(logistic.mcmcE1[,'L']))

logistic_o.postE1 <- exp(logistic.mcmcE1[,paste0('Lsubj[', 1:DATALIST_EXP1$S , ']')])
logistic_o.postE1 <- cbind(apply(logistic_o.postE1, 2, median), t(apply(logistic_o.postE1, 2, HDIofMCMC)))

dir.create('posteriors/')
save(list = c(ls(pattern = 'samplesE1'), ls(pattern = 'postE1')), file = 'posteriors/posteriorsExp1.RData')
