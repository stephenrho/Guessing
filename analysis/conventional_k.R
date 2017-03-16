### MORE CONVENTIONAL ESTIMATION OF K FOR EXPERIMENTS 1-4
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)

setwd("data/")

#### EXPERIMENT 1 -----
exp1 <- read.table('exp1/exp1_full')

exp1.hf <- ddply(exp1, c('ppt', 'SetSize'), summarize,
                 hr = mean(respChange[testChange == 1]),
                 far = mean(respChange[testChange == 0]))

pashler <- function(N, H, FA){
  FA = ifelse(FA < 1, FA, .99)
  a <- N*((H-F)/(1-F))
  ifelse(a < 0, 0, a)
}

exp1.hf$k <- with(exp1.hf, pashler(N = SetSize, H = hr, FA = far))

plot(aggregate(k ~ SetSize, data = exp1.hf, FUN = mean), ylim = c(0,6), type='b')

exp1.hf <- within(exp1.hf, {
  ppt <- factor(ppt)
  SetSize <- factor(SetSize)
})

aov1 <- aov(k ~ SetSize + Error(ppt/SetSize), data = exp1.hf)
summary(aov1)

aov1.2 <- aov(k ~ SetSize + Error(ppt/SetSize), data = subset(exp1.hf, SetSize != 2))
summary(aov1.2)

#### EXPERIMENT 2 -----
# given the use of choose functions in the formula for number of changes
# we use maximum likelihood to save time (rather than solving for k for each C)

# for this we make exp2 into an array with dim(4, 8, Ns) like Rouder et al. (2008)
# array[,,s] = subject's data
# array[1,,s] = data from 1 change condition (2 = 2 changes and so on...)
# array[p,1:4,s] = data from SS = 5 condition (5:8 = SS8)

exp2 <- read.table('exp2/exp2_full')

exp2agg <- ddply(exp2, .variables = c('ppt', 'SetSize', 'Cond', 'testChange'), summarise,
                 N = length(respChange),
                 Freq = sum(respChange))

Ns = unique(exp2agg$SetSize)
Ss = unique(exp2agg$ppt)
NCs = unique(exp2agg$Cond)
arrayExp2 <- array(NA, dim = c(length(NCs), length(Ns)*4, length(Ss)))

for (s in Ss){
  for (nc in NCs){
    frqs = c()
    for (n in Ns){
      fh = exp2agg$Freq[exp2agg$ppt == s & exp2agg$Cond == nc & exp2agg$SetSize == n]
      crm = exp2agg$N[exp2agg$ppt == s & exp2agg$Cond == nc & exp2agg$SetSize == n] - fh
      frqs = c(frqs, fh[2], crm[2], fh[1], crm[1])
    }
    arrayExp2[which(nc==NCs),,s] <- frqs
  }
}

## Estimate k via maximum likelihood
# multinomial -LL (used by model functions)
negLL <- function(y, p){
  ll = ifelse((y==0 & p==0) | p==0, 0, y*log(p))
  return(-sum(ll))
}

index=rep(1:length(Ns), each=4)

# model function
pashler_C <- function(k, u, N, C){
  # an extension of the model described by Pashler (1988, Percep. Psychophy.)
  # to apply across different numbers of possible changes
  ## Parameters ##
  # k = number of items in VWM
  # u = uninformed guessing rate
  # N = set size
  # C = number of possible changes 
  ### ### ### ### ### ### ### ### ### 
  k = ifelse(k > N - C, N, k)
  d = 1 - choose(N - k, C)/choose(N, C)
  
  p = 1:4
  p[1] = d + (1 - d)*u
  p[2] = 1 - p[1]
  p[3] = u
  p[4] = 1 - p[3]
  return(p)
}

# likelihood function (for our purposes only estimate separate k for each N)
ll.pashler_C <- function(par, y){
  # length(par) == 4 (k1:2, u1:2)
  ll=0
  for(j in 1:length(NCs)){ # for each number of changes
    for(i in 1:length(Ns)){ # for each set size
      p = pashler_C(k = par[i], u = par[length(Ns) + i], N = Ns[i], C = NCs[j])
      ll = ll + negLL(y[j, index==i], p)
    }
  }
  if(any(c(par < rep(0,4), par > c(rep(max(Ns), 2), 1, 1)))){
    ll = ll + 10000 # penalty for going out of range
  }
  return(ll)
}

nPar = 4
nRep = 100

out <- matrix(nrow = length(Ss), ncol = nPar+1)
for (s in 1:length(Ss)){ 
  # for each participant
  # start the routine from different starting points and only keep lowest ll
  current.ll = Inf # initialize log lik counter (stores the current candidate)
  for (rep in 1:nRep){
    par = runif(n = nPar, min = c(0,0,0,0), max = c(8,8,1,1))
    res = optim(par, ll.pashler_C, y = arrayExp2[,,s])
    if (res$value < current.ll){ 
      # if this is the smallest ll so far
      current.ll <- res$value
      out[s,] <- c(res$par, res$value)
    }
  }
}

apply(X = out, 2, mean)

t.test(out[,1], out[,2], paired = T)

# double check with single change block where conventional pashler formula can be applied
exp2.hf <- ddply(subset(exp2, Cond == 1), c('ppt', 'SetSize'), summarize,
                 hr = mean(respChange[testChange == 1]),
                 far = mean(respChange[testChange == 0]))

exp2.hf$k <- with(exp2.hf, pashler(N = SetSize, H = hr, FA = far))

plot(aggregate(k ~ SetSize, data = exp2.hf, FUN = mean), ylim = c(0,6), type='b')

exp2.hf <- within(exp2.hf, {
  ppt <- factor(ppt)
  SetSize <- factor(SetSize)
})

aov2 <- aov(k ~ SetSize + Error(ppt/SetSize), data = exp2.hf)
summary(aov2)

#### EXPERIMENT 3 -----
exp3 <- read.table('exp3/exp3_full')

exp3.hf <- ddply(exp3, c('ppt', 'SetSize'), summarize,
                 hr = mean(respChange[testChange == 1]),
                 far = mean(respChange[testChange == 0]))

rev.pashler <- function(N, H, FA){
  H = ifelse(H > .01, H, 0.01)
  a <- N*((H - FA)/ H)
  ifelse(a < 0, 0, a)
}

exp3.hf$k <- with(exp3.hf, rev.pashler(N = SetSize, H = hr, FA = far))

plot(aggregate(k ~ SetSize, data = exp3.hf, FUN = mean), ylim = c(0,6), type='b')

exp3.hf <- within(exp3.hf, {
  ppt <- factor(ppt)
  SetSize <- factor(SetSize)
})

aov3 <- aov(k ~ SetSize + Error(ppt/SetSize), data = exp3.hf)
summary(aov3)

aov3.2 <- aov(k ~ SetSize + Error(ppt/SetSize), data = subset(exp3.hf, SetSize != 2))
summary(aov3.2)

#### EXPERIMENT 4 -----
exp4 <- read.table('exp4/exp4_full')

exp4.hf <- ddply(exp4, c('ppt', 'SetSize'), summarize,
                 hr = mean(respChange[testChange == 1]),
                 far = mean(respChange[testChange == 0]))

exp4.hf$k <- with(exp4.hf, rev.pashler(N = SetSize, H = hr, FA = far))

plot(aggregate(k ~ SetSize, data = exp4.hf, FUN = mean), ylim = c(0,6), type='b')

exp4.hf <- within(exp4.hf, {
  ppt <- factor(ppt)
  SetSize <- factor(SetSize)
})

aov4 <- aov(k ~ SetSize + Error(ppt/SetSize), data = exp4.hf)
summary(aov4)

aov4.2 <- aov(k ~ SetSize + Error(ppt/SetSize), data = subset(exp4.hf, SetSize != 2))
summary(aov4.2)

# in all cases there are significant main effects of set size on k

saveList = c(ls(pattern = '.hf'), ls(pattern = 'aov'), 'out')

save(list = saveList, file = '../conventional_k.RData')

