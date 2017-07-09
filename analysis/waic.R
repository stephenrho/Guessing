# calculate WAIC for each experiment and model fit
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

load('guessingDataLists.RData')

source('waicFunctions.R')

## EXPERIMENT 1 ----

waic_dataME1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(dataM.samplesE1)), likFun = lik_dataM, predFunc = NULL)

waic_infE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(inf.samplesE1)), likFun = lik_guessing, predFunc = informed_2)

waic_uninfE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(uninf.samplesE1)), likFun = lik_guessing, predFunc = uninformed_2)

waic_inf_vkE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(inf_vk.samplesE1)), likFun = lik_guessing, predFunc = informed_2)

waic_uninf_vkE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(uninf_vk.samplesE1)), likFun = lik_guessing, predFunc = uninformed_2)

waic_mixtureE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(mixture.samplesE1)), likFun = lik_mixture, predFunc = mixture_2)

waic_logisticE1 = getWAIC(dataList = DATALIST_EXP1, samples = as.matrix(as.mcmc(logistic.samplesE1)), likFun = lik_logistic, predFunc = logistic_2)

## EXPERIMENT 2 ----

waic_dataME2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(dataM.samplesE2)), likFun = lik_dataM, predFunc = NULL)

waic_infE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(inf.samplesE2)), likFun = lik_guessing, predFunc = informed_2)

waic_uninfE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(uninf.samplesE2)), likFun = lik_guessing, predFunc = uninformed_2)

waic_inf_vkE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(inf_vk.samplesE2)), likFun = lik_guessing, predFunc = informed_2)

waic_uninf_vkE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(uninf_vk.samplesE2)), likFun = lik_guessing, predFunc = uninformed_2)

waic_mixtureE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(mixture.samplesE2)), likFun = lik_mixture, predFunc = mixture_2)

waic_logisticE2 = getWAIC(dataList = DATALIST_EXP2, samples = as.matrix(as.mcmc(logistic.samplesE2)), likFun = lik_logistic, predFunc = logistic_2)

## EXPERIMENT 3 ----

waic_dataME3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(dataM.samplesE3)), likFun = lik_dataM, predFunc = NULL)

waic_infE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(inf.samplesE3)), likFun = lik_guessing, predFunc = informed_3)

waic_uninfE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(uninf.samplesE3)), likFun = lik_guessing, predFunc = uninformed_3)

waic_inf_vkE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(inf_vk.samplesE3)), likFun = lik_guessing, predFunc = informed_3)

waic_uninf_vkE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(uninf_vk.samplesE3)), likFun = lik_guessing, predFunc = uninformed_3)

waic_mixtureE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(mixture.samplesE3)), likFun = lik_mixture, predFunc = mixture_3)

waic_logisticE3 = getWAIC(dataList = DATALIST_EXP3, samples = as.matrix(as.mcmc(logistic.samplesE3)), likFun = lik_logistic, predFunc = logistic_3)

## EXPERIMENT 4 ----

waic_dataME4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(dataM.samplesE4)), likFun = lik_dataM, predFunc = NULL)

waic_infE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(inf.samplesE4)), likFun = lik_guessing, predFunc = informed_3)

waic_uninfE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(uninf.samplesE4)), likFun = lik_guessing, predFunc = uninformed_3)

waic_inf_vkE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(inf_vk.samplesE4)), likFun = lik_guessing, predFunc = informed_3)

waic_uninf_vkE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(uninf_vk.samplesE4)), likFun = lik_guessing, predFunc = uninformed_3)

waic_mixtureE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(mixture.samplesE4)), likFun = lik_mixture, predFunc = mixture_3)

waic_logisticE4 = getWAIC(dataList = DATALIST_EXP4, samples = as.matrix(as.mcmc(logistic.samplesE4)), likFun = lik_logistic, predFunc = logistic_3)

# SAVE
save(list = ls(pattern = 'waic'), file = 'modelWAICs.RData')
