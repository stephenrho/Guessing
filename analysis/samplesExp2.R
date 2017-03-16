# run samples for experiment 2
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)
library(rjags)
library(R2jags)

# READ DATA AND MODELS ----------------------------------------------------
setwd('analysis/')

exp2 <- read.table(file = 'data/exp2/exp2_full')

# STRUCTURE DATA EXP 2 -----------------------------------------------------
DATALIST_EXP2 <- list(
  y = exp2$respChange,
  CHANGE = exp2$testChange,
  C = exp2$Cond,
  B_level = rep(1, nrow(exp2)),
  B_n = 1,
  N = exp2$SetSize,
  S = length(unique(exp2$ppt)),
  id = exp2$ppt,
  n = nrow(exp2),
  # for data model
  condLevel = exp2$Cond,
  N_c = length(unique(exp2$Cond)),
  ssLevel = as.numeric(as.factor(exp2$SetSize)),
  N_ss = length(unique(exp2$SetSize)),
  testLevel = as.numeric(as.factor(exp2$testChange))
)

# RUN CHAINS EXP 2 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD')

dataM.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(dataM.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD') 

inf.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/informedWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# variable k
inf_vk.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/informedWD_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf_vk.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# UNINFORMED
uninf.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/uninformedWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# variable k
uninf_vk.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/uninformedWD_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf_vk.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'P_og_subj')

mixture.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/mixtureWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(mixture.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# logistic RULE
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Lsubj')

logistic.samplesE2 <- jags.parallel(data = DATALIST_EXP2, inits = NULL, parameters.to.save = params, model.file = "models/logisticWD.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(logistic.samplesE2$BUGSoutput$summary)$Rhat < 1.1)

# EXTRACT POSTERIOR --------------------------------------------------------
source('guessingFunctions.R')

dataM.mcmcE2 <- as.matrix(as.mcmc(dataM.samplesE2))
inf.mcmcE2 <- as.matrix(as.mcmc(inf.samplesE2))
uninf.mcmcE2 <- as.matrix(as.mcmc(uninf.samplesE2))
mixture.mcmcE2 <- as.matrix(as.mcmc(mixture.samplesE2))
logistic.mcmcE2 <- as.matrix(as.mcmc(logistic.samplesE2))

dataM.postE2 <- logistic(dataM.mcmcE2[,3:18])

# extract estimated median f,h pairs and HDIs

fh.estimates.postE2 <- expand.grid(N = c(5,8), C = 1:4, 
                                   f.median=0, f.lower=0, f.upper=0,
                                   h.median=0, h.lower=0, h.upper=0)

for (nc in 1:4){
  for (n in 1:2){
    fcol <- paste0('R[1,', nc, ',', n, ']')
    hcol <- paste0('R[2,', nc, ',', n, ']')
    fh.estimates.postE2[fh.estimates.postE2$N == c(5,8)[n] & fh.estimates.postE2$C == c(1,2,3,4)[nc],3:8] <- c(median(dataM.postE2[,fcol]), HDIofMCMC(dataM.postE2[,fcol]), median(dataM.postE2[,hcol]), HDIofMCMC(dataM.postE2[,hcol]))
  }
}

# INFORMED
# extract posterior on natural scale
inf.postE2 <- cbind(inf.mcmcE2[,'K'],
                    logistic(inf.mcmcE2[,'A']),
                    logistic(inf.mcmcE2[,'U']))

# UNINFORMED
# extract posterior on natural scale
uninf.postE2 <- cbind(uninf.mcmcE2[,'K'],
                      logistic(uninf.mcmcE2[,'A']),
                      logistic(uninf.mcmcE2[,'U']))

# MIXTURE
mixture.postE2 <- cbind(mixture.mcmcE2[,'K'],
                        logistic(mixture.mcmcE2[,'A']),
                        logistic(mixture.mcmcE2[,'U']),
                        logistic(mixture.mcmcE2[,'POG']))

mixture_POG.postE2 <- logistic(mixture.mcmcE2[,paste0('P_og_subj[', 1:DATALIST_EXP2$S , ']')])
mixture_POG.postE2 <- cbind(apply(mixture_POG.postE2, 2, median), t(apply(mixture_POG.postE2, 2, HDIofMCMC)))

# logistic
logistic.postE2 <- cbind(logistic.mcmcE2[,'K'],
                      logistic(logistic.mcmcE2[,'A']),
                      logistic(logistic.mcmcE2[,'U']),
                      exp(logistic.mcmcE2[,'L']))

logistic_o.postE2 <- exp(logistic.mcmcE2[,paste0('Lsubj[', 1:DATALIST_EXP2$S , ']')])
logistic_o.postE2 <- cbind(apply(logistic_o.postE2, 2, median), t(apply(logistic_o.postE2, 2, HDIofMCMC)))

dir.create('posteriors/')
save(list = c(ls(pattern = 'samplesE2'), ls(pattern = 'postE2')), file = 'posteriors/posteriorsExp2.RData')
