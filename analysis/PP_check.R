# Posterior predictive check for
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())
library(plyr)
library(R2jags)

# READ EVERYTHING IN ----------------------------------------------------
# samples scripts should be run first to create posteriors
setwd('analysis/')

for (i in 1:4){
  load(paste0('posteriors/posteriorsExp', i, '.RData'))
}

rm(list=setdiff(ls(), ls(pattern = '^inf.')))

source('guessingFunctions.R')

#### EXPERIMENT 1 --------
# read data
exp1 <- read.table(file = 'data/exp1/exp1_full')
exp1$Cond <- as.factor(exp1$Cond)

exp1.agg <- ddply(exp1, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
              Ntrials = length(respChange),
              Nchange = sum(respChange))

mcmc.exp1 <- as.matrix(as.mcmc(inf.samplesE1))

# create matrix to hold pp samples
ppSamplesExp1 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), matrix(NA, nrow = 18, ncol = nrow(mcmc.exp1)+1))
ppSamplesExp1[,4] <- 60*abs(ppSamplesExp1[,2] + -1*as.numeric(!ppSamplesExp1[,3])) # n trials

for (i in 1:nrow(mcmc.exp1)){
  # sample parameters
  a = logistic(rnorm(1, mcmc.exp1[i, 'A'], mcmc.exp1[i, 'A_SD']))
  u_1 = logistic(rnorm(1, mcmc.exp1[i, 'U[1]'], mcmc.exp1[i, 'U_SD']))
  u_2 = logistic(rnorm(1, mcmc.exp1[i, 'U[2]'], mcmc.exp1[i, 'U_SD']))
  u_3 = logistic(rnorm(1, mcmc.exp1[i, 'U[3]'], mcmc.exp1[i, 'U_SD']))
  k = max(0, rnorm(1, mcmc.exp1[i, 'K'], mcmc.exp1[i, 'K_SD']))
  # sample n change responses
  p = genPredsInformedExp1(k, a, c(u_1, u_2, u_3), long = T)[,4]
  ppSamplesExp1[, 4+i] <- rbinom(n = 18, size = as.integer(ppSamplesExp1[,4]), prob = p)
}

# find central 95/50% of distribution
central95Exp1 <- t(apply(ppSamplesExp1[,5:ncol(ppSamplesExp1)], 1, function(x) HDIofMCMC(x, .95)))
central50Exp1 <- t(apply(ppSamplesExp1[,5:ncol(ppSamplesExp1)], 1, function(x) HDIofMCMC(x, .50)))

# what proportion of the data falls within the intervals?
exp1.agg$in95 = NA
exp1.agg$in50 = NA

for (s in unique(exp1.agg$ppt)){
  exp1.agg[exp1.agg$ppt == s,'in95'] <- as.integer(with(exp1.agg[exp1.agg$ppt == s,], Nchange >= central95Exp1[,1] & Nchange <= central95Exp1[,2]))
  exp1.agg[exp1.agg$ppt == s,'in50'] <- as.integer(with(exp1.agg[exp1.agg$ppt == s,], Nchange >= central50Exp1[,1] & Nchange <= central50Exp1[,2]))
}

#### EXPERIMENT 2 --------
# read data
exp2 <- read.table(file = 'data/exp2/exp2_full')
exp2$Cond <- as.factor(exp2$Cond)

exp2.agg <- ddply(exp2, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  Ntrials = length(respChange),
                  Nchange = sum(respChange))

mcmc.exp2 <- as.matrix(as.mcmc(inf.samplesE2))

# create matrix to hold pp samples
ppSamplesExp2 <- cbind(rep(c(5,8), each = 2*4), rep(1:4, each = 2), rep(c(0,1)), matrix(NA, nrow = 16, ncol = nrow(mcmc.exp2)+1))
ppSamplesExp2[,4] <- 34 # n trials

for (i in 1:nrow(mcmc.exp2)){
  # sample parameters
  a = logistic(rnorm(1, mcmc.exp2[i, 'A'], mcmc.exp2[i, 'A_SD']))
  u = logistic(rnorm(1, mcmc.exp2[i, 'U'], mcmc.exp2[i, 'U_SD']))
  k = max(0, rnorm(1, mcmc.exp2[i, 'K'], mcmc.exp2[i, 'K_SD']))
  # sample n change responses
  p = genPredsInformedExp2(k, a, u, long = T)[,4]
  ppSamplesExp2[, 4+i] <- rbinom(n = 16, size = as.integer(ppSamplesExp2[,4]), prob = p)
}

# find central 95/50% of distribution
central95Exp2 <- t(apply(ppSamplesExp2[,5:ncol(ppSamplesExp2)], 1, function(x) HDIofMCMC(x, .95)))
central50Exp2 <- t(apply(ppSamplesExp2[,5:ncol(ppSamplesExp2)], 1, function(x) HDIofMCMC(x, .50)))

# what proportion of the data falls within the intervals?
exp2.agg$in95 = NA
exp2.agg$in50 = NA

for (s in unique(exp2.agg$ppt)){
  exp2.agg[exp2.agg$ppt == s,'in95'] <- as.integer(with(exp2.agg[exp2.agg$ppt == s,], Nchange >= central95Exp2[,1] & Nchange <= central95Exp2[,2]))
  exp2.agg[exp2.agg$ppt == s,'in50'] <- as.integer(with(exp2.agg[exp2.agg$ppt == s,], Nchange >= central50Exp2[,1] & Nchange <= central50Exp2[,2]))
}

#### EXPERIMENT 3 --------
# read data
exp3 <- read.table(file = 'data/exp3/exp3_full')
exp3$Cond <- as.factor(exp3$Cond)

exp3.agg <- ddply(exp3, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  Ntrials = length(respChange),
                  Nchange = sum(respChange))

mcmc.exp3 <- as.matrix(as.mcmc(inf.samplesE3))

# create matrix to hold pp samples
ppSamplesExp3 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), matrix(NA, nrow = 18, ncol = nrow(mcmc.exp3)+1))
ppSamplesExp3[,4] <- 60*abs(ppSamplesExp3[,2] + -1*as.numeric(!ppSamplesExp3[,3])) # n trials

for (i in 1:nrow(mcmc.exp3)){
  # sample parameters
  a = logistic(rnorm(1, mcmc.exp3[i, 'A'], mcmc.exp3[i, 'A_SD']))
  u_1 = logistic(rnorm(1, mcmc.exp3[i, 'U[1]'], mcmc.exp3[i, 'U_SD']))
  u_2 = logistic(rnorm(1, mcmc.exp3[i, 'U[2]'], mcmc.exp3[i, 'U_SD']))
  u_3 = logistic(rnorm(1, mcmc.exp3[i, 'U[3]'], mcmc.exp3[i, 'U_SD']))
  k = max(0, rnorm(1, mcmc.exp3[i, 'K'], mcmc.exp3[i, 'K_SD']))
  # sample n change responses
  p = genPredsInformedExp3(k, a, c(u_1, u_2, u_3), long = T)[,4]
  ppSamplesExp3[, 4+i] <- rbinom(n = 18, size = as.integer(ppSamplesExp3[,4]), prob = p)
}

# find central 95/50% of distribution
central95Exp3 <- t(apply(ppSamplesExp3[,5:ncol(ppSamplesExp3)], 1, function(x) HDIofMCMC(x, .95)))
central50Exp3 <- t(apply(ppSamplesExp3[,5:ncol(ppSamplesExp3)], 1, function(x) HDIofMCMC(x, .50)))

# what proportion of the data falls within the intervals?
exp3.agg$in95 = NA
exp3.agg$in50 = NA

for (s in unique(exp3.agg$ppt)){
  exp3.agg[exp3.agg$ppt == s,'in95'] <- as.integer(with(exp3.agg[exp3.agg$ppt == s,], Nchange >= central95Exp3[,1] & Nchange <= central95Exp3[,2]))
  exp3.agg[exp3.agg$ppt == s,'in50'] <- as.integer(with(exp3.agg[exp3.agg$ppt == s,], Nchange >= central50Exp3[,1] & Nchange <= central50Exp3[,2]))
}

#### EXPERIMENT 4 --------
# read data
exp4 <- read.table(file = 'data/exp4/exp4_full')
exp4$Cond <- as.factor(exp4$Cond)

exp4.agg <- ddply(exp4, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  Ntrials = length(respChange),
                  Nchange = sum(respChange))

mcmc.exp4 <- as.matrix(as.mcmc(inf.samplesE4))

# create matrix to hold pp samples
ppSamplesExp4 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), matrix(NA, nrow = 18, ncol = nrow(mcmc.exp4)+1))
ppSamplesExp4[,4] <- 60*abs(ppSamplesExp4[,2] + -1*as.numeric(!ppSamplesExp4[,3])) # n trials

for (i in 1:nrow(mcmc.exp4)){
  # sample parameters
  a = logistic(rnorm(1, mcmc.exp4[i, 'A'], mcmc.exp4[i, 'A_SD']))
  u_1 = logistic(rnorm(1, mcmc.exp4[i, 'U[1]'], mcmc.exp4[i, 'U_SD']))
  u_2 = logistic(rnorm(1, mcmc.exp4[i, 'U[2]'], mcmc.exp4[i, 'U_SD']))
  u_3 = logistic(rnorm(1, mcmc.exp4[i, 'U[3]'], mcmc.exp4[i, 'U_SD']))
  k = max(0, rnorm(1, mcmc.exp4[i, 'K'], mcmc.exp4[i, 'K_SD']))
  # sample n change responses
  p = genPredsInformedExp3(k, a, c(u_1, u_2, u_3), long = T)[,4]
  ppSamplesExp4[, 4+i] <- rbinom(n = 18, size = as.integer(ppSamplesExp4[,4]), prob = p)
}

# find central 95/50% of distribution
central95Exp4 <- t(apply(ppSamplesExp4[,5:ncol(ppSamplesExp4)], 1, function(x) HDIofMCMC(x, .95)))
central50Exp4 <- t(apply(ppSamplesExp4[,5:ncol(ppSamplesExp4)], 1, function(x) HDIofMCMC(x, .50)))

# what proportion of the data falls within the intervals?
exp4.agg$in95 = NA
exp4.agg$in50 = NA

for (s in unique(exp4.agg$ppt)){
  exp4.agg[exp4.agg$ppt == s,'in95'] <- as.integer(with(exp4.agg[exp4.agg$ppt == s,], Nchange >= central95Exp4[,1] & Nchange <= central95Exp4[,2]))
  exp4.agg[exp4.agg$ppt == s,'in50'] <- as.integer(with(exp4.agg[exp4.agg$ppt == s,], Nchange >= central50Exp4[,1] & Nchange <= central50Exp4[,2]))
}

#### LOOK AT VALUES ----
mean(exp1.agg$in95)
mean(exp1.agg$in50)

mean(exp2.agg$in95)
mean(exp2.agg$in50)

mean(exp3.agg$in95)
mean(exp3.agg$in50)

mean(exp4.agg$in95)
mean(exp4.agg$in50)
