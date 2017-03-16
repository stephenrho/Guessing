# run samples for experiment 3
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)
library(rjags)
library(R2jags)

# READ DATA AND MODELS ----------------------------------------------------
setwd('analysis/')

exp3 <- read.table(file = 'data/exp3/exp3_full')
exp3$Cond <- as.factor(exp3$Cond)

# STRUCTURE DATA EXP 1 -----------------------------------------------------
DATALIST_EXP3 <- list(
  y = exp3$respChange,
  CHANGE = exp3$testChange,
  B_level = as.numeric(exp3$Cond),
  B_n = length(unique(exp3$Cond)),
  N = exp3$SetSize,
  S = length(unique(exp3$ppt)),
  id = exp3$ppt,
  n = nrow(exp3),
  # for data model
  condLevel = as.numeric(exp3$Cond),
  N_c = length(unique(exp3$Cond)),
  ssLevel = as.numeric(as.factor(exp3$SetSize)),
  N_ss = length(unique(exp3$SetSize)),
  testLevel = as.numeric(as.factor(exp3$testChange))
)

# RUN CHAINS EXP 3 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD')

dataM.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(dataM.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD') 

inf.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/informedSP.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

inf_vk.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/informedSP_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf_vk.samplesE3$BUGSoutput$summary)$Rhat < 1.1)
# K[1] is unstable. Increasing iterations will not solve

# UNINFORMED
uninf.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP.txt", n.chains = 5, n.iter = 55000*2, n.burnin = 5000, n.thin = 5*2)

all(as.data.frame(uninf.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

uninf_vk.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf_vk.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'P_og_subj')

mixture.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/mixtureSP.txt", n.chains = 5, n.iter = 55000, n.burnin = 5000, n.thin = 5)

all(as.data.frame(mixture.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

# logistic
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Lsubj')

logistic.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = 'models/logisticSP.txt', n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(logistic.samplesE3$BUGSoutput$summary)$Rhat < 1.1)

# EXTRACT POSTERIOR --------------------------------------------------------
source('guessingFunctions.R')

dataM.mcmcE3 <- as.matrix(as.mcmc(dataM.samplesE3))
inf.mcmcE3 <- as.matrix(as.mcmc(inf.samplesE3))
uninf.mcmcE3 <- as.matrix(as.mcmc(uninf.samplesE3))
mixture.mcmcE3 <- as.matrix(as.mcmc(mixture.samplesE3))
logistic.mcmcE3 <- as.matrix(as.mcmc(logistic.samplesE3))

dataM.postE3 <- logistic(dataM.mcmcE3[,!colnames(dataM.mcmcE3) %in% c('deviance', 'R_SD')])

# extract estimated median f,h pairs and HDIs
fh.estimates.postE3 <- expand.grid(N = c(2,5,8), P = c(0.3,0.5,0.7), 
                                   f.median=0, f.lower=0, f.upper=0,
                                   h.median=0, h.lower=0, h.upper=0)

for (p in 1:3){
  for (n in 1:3){
    fcol <- paste0('R[1,', p, ',', n, ']')
    hcol <- paste0('R[2,', p, ',', n, ']')
    fh.estimates.postE3[fh.estimates.postE3$N == c(2,5,8)[n] & fh.estimates.postE3$P == c(.3,.5,.7)[p],3:8] <- c(median(dataM.postE3[,fcol]), HDIofMCMC(dataM.postE3[,fcol]), median(dataM.postE3[,hcol]), HDIofMCMC(dataM.postE3[,hcol]))
  }
}

# INFORMED
# extract posterior on natural scale
inf.postE3 <- cbind(inf.mcmcE3[,'K'],
                    logistic(inf.mcmcE3[,'A']),
                    logistic(inf.mcmcE3[,'U[1]']),
                    logistic(inf.mcmcE3[,'U[2]']),
                    logistic(inf.mcmcE3[,'U[3]']))

# UNINFORMED
# extract posterior on natural scale
uninf.postE3 <- cbind(uninf.mcmcE3[,'K'],
                      logistic(uninf.mcmcE3[,'A']),
                      logistic(uninf.mcmcE3[,'U[1]']),
                      logistic(uninf.mcmcE3[,'U[2]']),
                      logistic(uninf.mcmcE3[,'U[3]']))

# MIXTURE
# extract posterior on natural scale
mixture.postE3 <- cbind(mixture.mcmcE3[,'K'],
                        logistic(mixture.mcmcE3[,'A']),
                        logistic(mixture.mcmcE3[,'U[1]']),
                        logistic(mixture.mcmcE3[,'U[2]']),
                        logistic(mixture.mcmcE3[,'U[3]']),
                        logistic(mixture.mcmcE3[,'POG']))

mixture_POG.postE3 <- logistic(mixture.mcmcE3[,paste0('P_og_subj[', 1:DATALIST_EXP3$S , ']')])
mixture_POG.postE3 <- cbind(apply(mixture_POG.postE3, 2, median), t(apply(mixture_POG.postE3, 2, HDIofMCMC)))

# logistic
logistic.postE3 <- cbind(logistic.mcmcE3[,'K'],
                      logistic(logistic.mcmcE3[,'A']),
                      logistic(logistic.mcmcE3[,'U[1]']),
                      logistic(logistic.mcmcE3[,'U[2]']),
                      logistic(logistic.mcmcE3[,'U[3]']),
                      exp(logistic.mcmcE3[,'L']))

logistic_o.postE3 <- exp(logistic.mcmcE3[,paste0('Lsubj[', 1:DATALIST_EXP3$S , ']')])
logistic_o.postE3 <- cbind(apply(logistic_o.postE3, 2, median), t(apply(logistic_o.postE3, 2, HDIofMCMC)))

dir.create('posteriors/')
save(list = c(ls(pattern = 'samplesE3'), ls(pattern = 'postE3')), file = 'posteriors/posteriorsExp3.RData')
