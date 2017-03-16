# run samples for experiment 4
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)
library(rjags)
library(R2jags)

# READ DATA AND MODELS ----------------------------------------------------
setwd('analysis/')

exp4 <- read.table(file = 'data/exp4/exp4_full')
exp4$Cond <- as.factor(exp4$Cond)

# STRUCTURE DATA EXP 1 -----------------------------------------------------
DATALIST_EXP4 <- list(
  y = exp4$respChange,
  CHANGE = exp4$testChange,
  B_level = as.numeric(exp4$Cond),
  B_n = length(unique(exp4$Cond)),
  N = exp4$SetSize,
  S = length(unique(exp4$ppt)),
  id = exp4$ppt,
  n = nrow(exp4),
  # for data model
  condLevel = as.numeric(exp4$Cond),
  N_c = length(unique(exp4$Cond)),
  ssLevel = as.numeric(as.factor(exp4$SetSize)),
  N_ss = length(unique(exp4$SetSize)),
  testLevel = as.numeric(as.factor(exp4$testChange))
)

# RUN CHAINS EXP 4 --------------------------------------------------------

# DATA MODEL
params <- c('R', 'R_SD')

dataM.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/dataModel.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(dataM.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

# INFORMED
params <- c('K', 'A', 'U', 'K_SD', 'A_SD', 'U_SD') 

inf.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/informedSP.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

inf_vk.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/informedSP_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(inf_vk.samplesE4$BUGSoutput$summary)$Rhat < 1.1)
# K[1] is unstable. Increasing iterations will not solve

# UNINFORMED
uninf.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

uninf_vk.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/uninformedSP_vk.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(uninf_vk.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

# MIXTURE
params <- c('K', 'A', 'U', 'POG', 'K_SD', 'A_SD', 'U_SD', 'POG_SD', 'P_og_subj')

mixture.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/mixtureSP.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(mixture.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

# logistic RULE
params <- c('K', 'A', 'U', 'L', 'K_SD', 'A_SD', 'U_SD', 'L_SD', 'Lsubj')

logistic.samplesE4 <- jags.parallel(data = DATALIST_EXP4, inits = NULL, parameters.to.save = params, model.file = "models/logisticSP.txt", n.chains = 5, n.iter = 15000, n.burnin = 5000, n.thin = 1)

all(as.data.frame(logistic.samplesE4$BUGSoutput$summary)$Rhat < 1.1)

# EXTRACT POSTERIOR --------------------------------------------------------
source('guessingFunctions.R')

dataM.mcmcE4 <- as.matrix(as.mcmc(dataM.samplesE4))
inf.mcmcE4 <- as.matrix(as.mcmc(inf.samplesE4))
uninf.mcmcE4 <- as.matrix(as.mcmc(uninf.samplesE4))
mixture.mcmcE4 <- as.matrix(as.mcmc(mixture.samplesE4))
logistic.mcmcE4 <- as.matrix(as.mcmc(logistic.samplesE4))

dataM.postE4 <- logistic(dataM.mcmcE4[,!colnames(dataM.mcmcE4) %in% c('deviance', 'R_SD')])

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

dir.create('posteriors/')
save(list = c(ls(pattern = 'samplesE4'), ls(pattern = 'postE4')), file = 'posteriors/posteriorsExp4.RData')
