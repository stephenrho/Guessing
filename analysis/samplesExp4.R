# run samples for experiment 4
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

# RUN CHAINS EXP 4 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD', 'Rsubj')

dataM.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

dataM.summaryE4 = as.data.frame(dataM.samplesE4$BUGSoutput$summary)
dataM.summaryE4[dataM.summaryE4$Rhat>1.1,]

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD', 'Ksubj', 'Asubj', 'Usubj') 

inf.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/informedSP.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

inf.summaryE4 = as.data.frame(inf.samplesE4$BUGSoutput$summary)
inf.summaryE4[inf.summaryE4$Rhat>1.1,]

DATALIST_EXP4['K1_sd'] = .1 # without prior constraint K[1] is unstable
inf_vk.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/informedSP_vk.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

inf_vk.summaryE4 = as.data.frame(inf_vk.samplesE4$BUGSoutput$summary)
inf_vk.summaryE4[inf_vk.summaryE4$Rhat>1.1,]

# UNINFORMED
uninf.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP.txt", n.chains = 4, n.iter = 105000, n.burnin = 5000, n.thin = 20)

uninf.summaryE4 = as.data.frame(uninf.samplesE4$BUGSoutput$summary)
uninf.summaryE4[uninf.summaryE4$Rhat>1.1,]

uninf_vk.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vk.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

uninf_vk.summaryE4 = as.data.frame(uninf_vk.samplesE4$BUGSoutput$summary)
uninf_vk.summaryE4[uninf_vk.summaryE4$Rhat>1.1,]

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'Ksubj', 'Asubj', 'Usubj', 'P_og_subj')

mixture.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/mixtureSP.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

mixture.summaryE4 = as.data.frame(mixture.samplesE4$BUGSoutput$summary)
mixture.summaryE4[mixture.summaryE4$Rhat>1.1,]

# logistic RULE
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Ksubj', 'Asubj', 'Usubj', 'Lsubj')

logistic.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = k_initFunc1, parameters.to.save = params, model.file = "models/logisticSP.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 1)

logistic.summaryE4 = as.data.frame(logistic.samplesE4$BUGSoutput$summary)
logistic.summaryE4[logistic.summaryE4$Rhat>1.1,]

# Allow u to vary by set size and condition
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD')

uninf_vu.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vu.txt", n.chains = 4, n.iter = 105000, n.burnin = 5000, n.thin = 20)

uninf_vu.summaryE4 = as.data.frame(uninf_vu.samplesE4$BUGSoutput$summary)
uninf_vu.summaryE4[uninf_vu.summaryE4$Rhat>1.1,]

# EXTRACT POSTERIOR --------------------------------------------------------

dataM.mcmcE4 <- as.matrix(as.mcmc(dataM.samplesE4))
inf.mcmcE4 <- as.matrix(as.mcmc(inf.samplesE4))
uninf.mcmcE4 <- as.matrix(as.mcmc(uninf.samplesE4))
mixture.mcmcE4 <- as.matrix(as.mcmc(mixture.samplesE4))
logistic.mcmcE4 <- as.matrix(as.mcmc(logistic.samplesE4))

dataM.postE4 <- logistic(dataM.mcmcE4[,3:20])

# extract estimated median f,h pairs and HDIs
fh.estimates.postE4 <- expand.grid(N = c(2,5,8), P = c(0.3,0.5,0.7), 
                                   f.median=0, f.lower=0, f.upper=0,
                                   h.median=0, h.lower=0, h.upper=0)

for (p in 1:3){
  for (n in 1:3){
    fcol <- paste0('R[1,', p, ',', n, ']')
    hcol <- paste0('R[2,', p, ',', n, ']')
    fh.estimates.postE4[fh.estimates.postE4$N == c(2,5,8)[n] & fh.estimates.postE4$P == c(.3,.5,.7)[p],3:8] <- c(median(dataM.postE4[,fcol]), HDIofMCMC(dataM.postE4[,fcol]), median(dataM.postE4[,hcol]), HDIofMCMC(dataM.postE4[,hcol]))
  }
}

# INFORMED
# extract posterior on natural scale
inf.postE4 <- cbind(inf.mcmcE4[,'K'],
                    logistic(inf.mcmcE4[,'A']),
                    logistic(inf.mcmcE4[,'U[1]']),
                    logistic(inf.mcmcE4[,'U[2]']),
                    logistic(inf.mcmcE4[,'U[3]']))

# UNINFORMED
# extract posterior on natural scale
uninf.postE4 <- cbind(uninf.mcmcE4[,'K'],
                      logistic(uninf.mcmcE4[,'A']),
                      logistic(uninf.mcmcE4[,'U[1]']),
                      logistic(uninf.mcmcE4[,'U[2]']),
                      logistic(uninf.mcmcE4[,'U[3]']))

# MIXTURE
# extract posterior on natural scale
mixture.postE4 <- cbind(mixture.mcmcE4[,'K'],
                        logistic(mixture.mcmcE4[,'A']),
                        logistic(mixture.mcmcE4[,'U[1]']),
                        logistic(mixture.mcmcE4[,'U[2]']),
                        logistic(mixture.mcmcE4[,'U[3]']),
                        logistic(mixture.mcmcE4[,'POG']))

mixture_POG.postE4 <- logistic(mixture.mcmcE4[,paste0('P_og_subj[', 1:DATALIST_EXP4$S , ']')])
mixture_POG.postE4 <- cbind(apply(mixture_POG.postE4, 2, median), t(apply(mixture_POG.postE4, 2, HDIofMCMC)))

# logistic
logistic.postE4 <- cbind(logistic.mcmcE4[,'K'],
                      logistic(logistic.mcmcE4[,'A']),
                      logistic(logistic.mcmcE4[,'U[1]']),
                      logistic(logistic.mcmcE4[,'U[2]']),
                      logistic(logistic.mcmcE4[,'U[3]']),
                      exp(logistic.mcmcE4[,'L']))

logistic_o.postE4 <- exp(logistic.mcmcE4[,paste0('Lsubj[', 1:DATALIST_EXP4$S , ']')])
logistic_o.postE4 <- cbind(apply(logistic_o.postE4, 2, median), t(apply(logistic_o.postE4, 2, HDIofMCMC)))

save(list = c(ls(pattern = 'samplesE4'), ls(pattern = 'postE4')), file = 'posteriors/posteriorsExp4.RData')

