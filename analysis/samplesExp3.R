# run samples for experiment 3
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)
library(rjags)
library(R2jags)

# READ DATA AND FUNCTIONS ----------------------------------------------------
setwd('analysis/')

load('guessingDataLists.RData')

source('guessingFunctions.R')

# RUN CHAINS EXP 3 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD', 'Rsubj')

dataM.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

dataM.summaryE3 = as.data.frame(dataM.samplesE3$BUGSoutput$summary)
dataM.summaryE3[dataM.summaryE3$Rhat>1.1,]

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD', 'Ksubj', 'Asubj', 'Usubj') 

inf.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/informedSP.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

inf.summaryE3 = as.data.frame(inf.samplesE3$BUGSoutput$summary)
inf.summaryE3[inf.summaryE3$Rhat>1.1,]

DATALIST_EXP3['K1_sd'] = .1 # without prior constraint K[1] is unstable
inf_vk.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/informedSP_vk.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

inf_vk.summaryE3 = as.data.frame(inf_vk.samplesE3$BUGSoutput$summary)
inf_vk.summaryE3[inf_vk.summaryE3$Rhat>1.1,]

# UNINFORMED
uninf.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP.txt", n.chains = 4, n.iter = 145000, n.burnin = 5000, n.thin = 30)

uninf.summaryE3 = as.data.frame(uninf.samplesE3$BUGSoutput$summary)
uninf.summaryE3[uninf.summaryE3$Rhat>1.1,]

uninf_vk.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vk.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

uninf_vk.summaryE3 = as.data.frame(uninf_vk.samplesE3$BUGSoutput$summary)
uninf_vk.summaryE3[uninf_vk.summaryE3$Rhat>1.1,]

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'Ksubj', 'Asubj', 'Usubj', 'P_og_subj')

mixture.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/mixtureSP.txt", n.chains = 4, n.iter = 30000, n.burnin = 5000, n.thin = 5)

mixture.summaryE3 = as.data.frame(mixture.samplesE3$BUGSoutput$summary)
mixture.summaryE3[mixture.summaryE3$Rhat>1.1,]

# logistic
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Ksubj', 'Asubj', 'Usubj', 'Lsubj')

logistic.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = k_initFunc1, parameters.to.save = params, model.file = 'models/logisticSP.txt', n.chains = 4, n.iter = 105000, n.burnin = 5000, n.thin = 20)

logistic.summaryE3 = as.data.frame(logistic.samplesE3$BUGSoutput$summary)
logistic.summaryE3[logistic.summaryE3$Rhat>1.1,]

# Allow u to vary by set size and condition
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD')

uninf_vu.samplesE3 <- jags.parallel(data = DATALIST_EXP3, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vu.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

uninf_vu.summaryE3 = as.data.frame(uninf_vu.samplesE3$BUGSoutput$summary)
uninf_vu.summaryE3[uninf_vu.summaryE3$Rhat>1.1,]

# EXTRACT POSTERIOR --------------------------------------------------------

dataM.mcmcE3 <- as.matrix(as.mcmc(dataM.samplesE3))
inf.mcmcE3 <- as.matrix(as.mcmc(inf.samplesE3))
uninf.mcmcE3 <- as.matrix(as.mcmc(uninf.samplesE3))
mixture.mcmcE3 <- as.matrix(as.mcmc(mixture.samplesE3))
logistic.mcmcE3 <- as.matrix(as.mcmc(logistic.samplesE3))

dataM.postE3 <- logistic(dataM.mcmcE3[,3:20])

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

save(list = c(ls(pattern = 'samplesE3'), ls(pattern = 'postE3')), file = 'posteriors/posteriorsExp3.RData')

