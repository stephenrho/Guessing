#### LOOK AT PARTICIPANT LEVEL ESTIMATES OF K
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())
library(plyr)
library(R2jags)

# READ EVERYTHING IN ----------------------------------------------------
setwd('analysis/')

for (i in 1:4){
  load(paste0('posteriors/posteriorsExp', i, '.RData'))
}

source('guessingFunctions.R')

# to extract median and hdis of individual k samples
getIndK = function(samples){
  samples = samples[,grep(pattern = 'Ksubj', x = colnames(samples))]
  
  m = apply(samples, 2, median)
  hdi = t(apply(samples, 2, HDIofMCMC))
  tmp = cbind(m, hdi)
  colnames(tmp) = c('m', 'l', 'u')
  return(as.data.frame(tmp))
}

# 1
infE1 = getIndK(as.matrix(as.mcmc(inf.samplesE1)))
any(infE1$l < 5 & infE1$u > 5)

uninfE1 = getIndK(as.matrix(as.mcmc(uninf.samplesE1)))
any(uninfE1$l < 5 & uninfE1$u > 5)

# 2
infE2 = getIndK(as.matrix(as.mcmc(inf.samplesE2)))
any(infE2$l < 5 & infE2$u > 5)

uninfE2 = getIndK(as.matrix(as.mcmc(uninf.samplesE2)))
any(uninfE2$l < 5 & uninfE2$u > 5)

uninfE2[which(uninfE2$l < 5 & uninfE2$u > 5),]
#           m        l        u
# Ksubj[21] 5.003087 3.617699 6.095437
# Ksubj[27] 4.167613 3.999598 5.263729

# 3
infE3 = getIndK(as.matrix(as.mcmc(inf.samplesE3)))
any(infE3$l < 5 & infE3$u > 5)

infE3[which(infE3$l < 5 & infE3$u > 5),]
#          m        l        u
# Ksubj[4] 4.699376 4.306464 5.039887

uninfE3 = getIndK(as.matrix(as.mcmc(uninf.samplesE3)))
any(uninfE3$l < 5 & uninfE3$u > 5)

uninfE3[which(uninfE3$l < 5 & uninfE3$u > 5),]
#          m      l        u
# Ksubj[4] 5.129162 4.5195 5.643892
# Ksubj[8] 4.215467 3.6391 5.197484

# 4
infE4 = getIndK(as.matrix(as.mcmc(inf.samplesE4)))
any(infE4$l < 5 & infE4$u > 5)

infE4[which(infE4$l < 5 & infE4$u > 5),]
#           m        l        u
# Ksubj[12] 4.664794 4.238182 5.473598

uninfE4 = getIndK(as.matrix(as.mcmc(uninf.samplesE4)))
any(uninfE4$l < 5 & uninfE4$u > 5)

uninfE4[which(uninfE4$l < 5 & uninfE4$u > 5),]
#           m        l        u
# Ksubj[11] 4.517973 3.643459 5.836098
# Ksubj[12] 4.857703 4.349461 5.976617
