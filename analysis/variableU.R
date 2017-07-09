# Look at estimated guessing rates from the unconstrained U model suggested by reviewer
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())
library(plyr)
library(R2jags)
library(vioplot)

# READ EVERYTHING IN ----------------------------------------------------
setwd('analysis/')

for (i in 1:4){
  load(paste0('posteriors/posteriorsExp', i, '.RData'))
}

# settings for figures
dir.create('vary_u_figs/')
hw = c(5,4)

source('guessingFunctions.R')

# U[ss, cond]

uninf_vu.postE1 = as.matrix(as.mcmc(uninf_vu.samplesE1))
uninf_vu.postE1 = plogis(uninf_vu.postE1[,grep('U\\[', colnames(uninf_vu.postE1))])

uninf_vu.postE2 = as.matrix(as.mcmc(uninf_vu.samplesE2))
uninf_vu.postE2 = plogis(uninf_vu.postE2[,grep('U\\[', colnames(uninf_vu.postE2))])

uninf_vu.postE3 = as.matrix(as.mcmc(uninf_vu.samplesE3))
uninf_vu.postE3 = plogis(uninf_vu.postE3[,grep('U\\[', colnames(uninf_vu.postE3))])

uninf_vu.postE4 = as.matrix(as.mcmc(uninf_vu.samplesE4))
uninf_vu.postE4 = plogis(uninf_vu.postE4[,grep('U\\[', colnames(uninf_vu.postE4))])

# PLOT ESTIMATES OF U BY SET SIZE ------
# exp 1
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    vioplot(uninf_vu.postE1[,paste0('U[', s, ',', u, ']')], at=s, col='violetred4', add=T)
  }
  axis(2)
  text(x = 3.2, y = .8, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 1 - variable u', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('u', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp1_vu.pdf', height = hw[1], width = hw[2])
dev.off()

# exp 2
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=rep(0,4))

plot(NA, xlim = c(.7, 2.3), ylim=c(0,1), xlab='', ylab='', axes=F)
for (s in 1:3){
  vioplot(uninf_vu.postE2[,paste0('U[', s, ',', 1, ']')], at=s, col='turquoise4', add=T)
}
axis(2)
axis(1, at = 1:2, labels = c(5,8))
mtext('Experiment 2 - variable u', 3)
mtext('N', 1, line = 2.5)
mtext('u', 2, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp2_vu.pdf', height = hw[1], width = hw[2])
dev.off()

# Exp 3
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    vioplot(uninf_vu.postE3[,paste0('U[', s, ',', u, ']')], at=s, col='tomato4', add=T)
  }
  axis(2)
  text(x = 3.2, y = .2, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 3 - variable u', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('u', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp3_vu.pdf', height = hw[1], width = hw[2])
dev.off()

# Exp 4
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    vioplot(uninf_vu.postE4[,paste0('U[', s, ',', u, ']')], at=s, col='steelblue4', add=T)
  }
  axis(2)
  text(x = 3.2, y = .2, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 4 - variable u', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('u', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp4_vu.pdf', height = hw[1], width = hw[2])
dev.off()

### COMPARE TO ESTIMATES OF G FROM INFORMED MODEL CHAINS -----

calcG_wd = function(k, u, N, C){
  d = ifelse(k >= N - C + 1, 1, 1 - choose(N - k, C)/choose(N, C))
  g = ((1 - d)*u)/((1 - d)*u + (1 - u))
  return(g)
}

calcG_sp = function(k, u, N){
  d = min(k/N, 1)
  g = u/(u + (1 - d)*(1 - u))
  return(g)
}

# create matrices of g samples
g.postE1 = c()
for (n in c(2,5,8)){
  for (p in 3:5){
    g.postE1 = cbind(g.postE1, calcG_wd(k = inf.postE1[,1], u = inf.postE1[,p], N = n, C=1))
  }
}
colnames(g.postE1) = c("G[1,1]", "G[1,2]", "G[1,3]", "G[2,1]", "G[2,2]", "G[2,3]", "G[3,1]", "G[3,2]", "G[3,3]")

g.postE2 = c()
for (n in c(5,8)){
  g.postE2 = cbind(g.postE2, calcG_wd(k = inf.postE2[,1], u = inf.postE2[,3], N = n, C=1))
}
colnames(g.postE2) = c("G[1,1]", "G[2,1]")

g.postE3 = c()
for (n in c(2,5,8)){
  for (p in 3:5){
    g.postE3 = cbind(g.postE3, calcG_sp(k = inf.postE3[,1], u = inf.postE3[,p], N = n))
  }
}
colnames(g.postE3) = c("G[1,1]", "G[1,2]", "G[1,3]", "G[2,1]", "G[2,2]", "G[2,3]", "G[3,1]", "G[3,2]", "G[3,3]")

g.postE4 = c()
for (n in c(2,5,8)){
  for (p in 3:5){
    g.postE4 = cbind(g.postE4, calcG_sp(k = inf.postE4[,1], u = inf.postE4[,p], N = n))
  }
}
colnames(g.postE4) = c("G[1,1]", "G[1,2]", "G[1,3]", "G[2,1]", "G[2,2]", "G[2,3]", "G[3,1]", "G[3,2]", "G[3,3]")


# plot g estimates from informed model
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    if(s == 1){
      points(x=s, y=0)
    } else{
      vioplot(g.postE1[,paste0('G[', s, ',', u, ']')], at=s, col='violet', add=T)
    }
  }
  axis(2)
  text(x = 3.2, y = .8, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 1 - informed g', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('g', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp1_g.pdf', height = hw[1], width = hw[2])
dev.off()

# exp 2
par(mfrow=c(1,1), mar=c(4,4,2,2), oma=rep(0,4))

plot(NA, xlim = c(.7, 2.3), ylim=c(0,1), xlab='', ylab='', axes=F)
for (s in 1:3){
  vioplot(g.postE2[,paste0('G[', s, ',', 1, ']')], at=s, col='turquoise', add=T)
}
axis(2)
axis(1, at = 1:2, labels = c(5,8))
mtext('Experiment 2 - informed g', 3)
mtext('N', 1, line = 2.5)
mtext('g', 2, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp2_g.pdf', height = hw[1], width = hw[2])
dev.off()

# Exp 3
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    if(s == 1){
      points(x=s, y=1)
    } else{
      vioplot(g.postE3[,paste0('G[', s, ',', u, ']')], at=s, col='tomato', add=T)
    }
  }
  axis(2)
  text(x = 3.2, y = .2, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 3 - informed g', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('g', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp3_g.pdf', height = hw[1], width = hw[2])
dev.off()

# Exp 4
par(mfrow=c(3,1), mar=c(.6,.0,0,0), oma=c(4,4,2,2))

for (u in 1:3){
  plot(NA, xlim = c(.7, 3.3), ylim=c(0,1), xlab='', ylab='', axes=F)
  for (s in 1:3){
    if(s == 1){
      points(x=s, y=1)
    } else{
      vioplot(g.postE4[,paste0('G[', s, ',', u, ']')], at=s, col='steelblue', add=T)
    }
  }
  axis(2)
  text(x = 3.2, y = .2, labels = paste0('p(change) = ', c('0.3', '0.5', '0.7')[u]), adj = 1)
}
axis(1, at = 1:3, labels = c(2,5,8))
mtext('Experiment 4 - informed g', 3, outer = T)
mtext('N', 1, outer = T, line = 2.5)
mtext('g', 2, outer = T, line = 2.5)

dev.copy(pdf, 'vary_u_figs/exp4_g.pdf', height = hw[1], width = hw[2])
dev.off()
