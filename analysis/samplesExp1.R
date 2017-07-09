# run samples for experiment 1
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

# RUN CHAINS EXP 1 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD', 'Rsubj')

dataM.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

dataM.summaryE1 = as.data.frame(dataM.samplesE1$BUGSoutput$summary)
dataM.summaryE1[dataM.summaryE1$Rhat>1.1,]

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD', 'Ksubj', 'Asubj', 'Usubj') 

inf.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc1, parameters.to.save = params, model.file = "models/informedWD.txt", n.chains = 4, n.iter = 55000, n.burnin = 5000, n.thin = 10)

inf.summaryE1 = as.data.frame(inf.samplesE1$BUGSoutput$summary)
inf.summaryE1[inf.summaryE1$Rhat>1.1,]

# variable k
DATALIST_EXP1['K1_sd'] = .1 # without prior constraint K[1] is unstable
inf_vk.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc3, parameters.to.save = params, model.file = "models/informedWD_vk.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

inf_vk.summaryE1 = as.data.frame(inf_vk.samplesE1$BUGSoutput$summary)
inf_vk.summaryE1[inf_vk.summaryE1$Rhat>1.1,]

# UNINFORMED
uninf.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc1, parameters.to.save = params, model.file = "models/uninformedWD.txt", n.chains = 4, n.iter = 105000, n.burnin = 5000, n.thin = 20)
# additional iterations and thinning run to converge

uninf.summaryE1 = as.data.frame(uninf.samplesE1$BUGSoutput$summary)
uninf.summaryE1[uninf.summaryE1$Rhat>1.1,]

# variable k
uninf_vk.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc3, parameters.to.save = params, model.file = "models/uninformedWD_vk2.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

uninf_vk.summaryE1 = as.data.frame(uninf_vk.samplesE1$BUGSoutput$summary)
uninf_vk.summaryE1[uninf_vk.summaryE1$Rhat>1.1,]

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'Ksubj', 'Asubj', 'Usubj', 'P_og_subj')

mixture.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc1, parameters.to.save = params, model.file = "models/mixtureWD.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

mixture.summaryE1 = as.data.frame(mixture.samplesE1$BUGSoutput$summary)
mixture.summaryE1[mixture.summaryE1$Rhat>1.1,]

# logistic RULE
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Ksubj', 'Asubj', 'Usubj', 'Lsubj')

logistic.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc1, parameters.to.save = params, model.file = "models/logisticWD.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

logistic.summaryE1 = as.data.frame(logistic.samplesE1$BUGSoutput$summary)
logistic.summaryE1[logistic.summaryE1$Rhat>1.1,]

# Allow u to vary by set size and condition
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD')

uninf_vu.samplesE1 <- jags.parallel(data = DATALIST_EXP1, inits = k_initFunc1, parameters.to.save = params, model.file = "models/uninformedWD_vu.txt", n.chains = 4, n.iter = 105000, n.burnin = 5000, n.thin = 20)

uninf_vu.summaryE1 = as.data.frame(uninf_vu.samplesE1$BUGSoutput$summary)
uninf_vu.summaryE1[uninf_vu.summaryE1$Rhat>1.1,]

# EXTRACT POSTERIOR --------------------------------------------------------

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

save(list = c(ls(pattern = 'samplesE1'), ls(pattern = 'postE1')), file = 'posteriors/posteriorsExp1.RData')
