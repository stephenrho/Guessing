# luminance and p(response = change)
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

# follows up on a reviewer comment to assess the role of the difference in average luminance
# between study and test arrays in determining the probability of respoding change

rm(list=ls())

library(plyr)
#library(rjags)
#library(R2jags)
library(rstanarm)
library(loo)

# COLORS
# psychopy colors given from -1 to 1 (0 = grey; RGB=255/2)
colorsSetPsychopy = rbind(
  c(-1,-1,-1), # 'black'
  c( 1, 1, 1), # 'white'
  c( 1,-1,-1), # 'red'
  c(-1,-1, 1), # 'blue', 
  c(-1, 0,-1), # 'green', 
  c(1, 1,-1), # 'yellow', 
  c(1,.3,-1), # 'orange', 
  c(-1, 1, 1), # 'cyan', 
  c(0,-1, 0), # 'purple', 
  c(-1, 0, 0), #'darkbluegreen' / teal...
  # experiment 2 only
  c(1,-1, 1), # fuchsia
  c(-1, 1,-1), # lime
  c(0,-1,-1) # maroon
)

# scale to 0-1
colorsSet = (colorsSetPsychopy+1)/2

# https://en.wikipedia.org/wiki/Relative_luminance
lumScale = c(.2126, .7152, .0722)

relativeLumSet = colorsSet %*% lumScale

extractMeanRelLum = function(colorVec, N, split = ''){
  # covert the number strings to color values to mean relative luminance
  out = rep(0, length(colorVec))
  for (i in 1:length(colorVec)){
    vec = as.numeric(unlist(strsplit(as.character(colorVec[i]), split = split)))
    if (length(vec) != N[i]){ # problems with 0 as first item
      vec = c(0, vec)
    }
    
    out[i] = mean(relativeLumSet[vec+1])
  }
  return(out)
}

# READ DATA ----------------------------------------------------
try(setwd('~/Dropbox/Guessing/write up/JAGS_analysis/'), silent = T)
try(setwd('C:/Users/rhodessp/Dropbox/Guessing/write up/JAGS_analysis/'), silent = T)

## Experiment 1
exp1 <- read.table(file = 'data/exp1/exp1_full')
exp1 = within(exp1, {
  Cond = as.factor(Cond)
  SetSize = as.factor(SetSize)
  ppt = as.factor(ppt)
})

exp1$avStudyLum = extractMeanRelLum(colorVec = exp1$study, N = exp1$SetSize)
exp1$avTestLum = extractMeanRelLum(colorVec = exp1$test, N = exp1$SetSize)

exp1$avLumDiff = with(exp1, avTestLum - avStudyLum)

exp1$absLumDiff = abs(exp1$avLumDiff)

options(contrasts = c('contr.sum', 'contr.sum'))

exp1_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp1, testChange==1), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp1_m1_loo=loo(exp1_m1)
exp1_m1_waic=waic(exp1_m1)

# include luminance as a predictor
exp1_m2 = update(exp1_m1, .~. + absLumDiff)

exp1_m2_loo=loo(exp1_m2)
exp1_m2_waic=waic(exp1_m2)

compare(exp1_m1_loo, exp1_m2_loo)
compare(exp1_m1_waic, exp1_m2_waic)

# WAIC is lower for the model excluding luminance difference as a predictor
summary(exp1_m2) # p(resp=change) is not really moderated by luminance difference

## Experiment 2
exp2 <- read.table(file = 'data/exp2/exp2_full')
exp2 = within(exp2, {
  Cond = as.factor(Cond)
  SetSize = as.factor(SetSize)
  ppt = as.factor(ppt)
})

exp2$avStudyLum = extractMeanRelLum(colorVec = exp2$study, N = exp2$SetSize, split = '-')
exp2$avTestLum = extractMeanRelLum(colorVec = exp2$test, N = exp2$SetSize, split = '-')

exp2$avLumDiff = with(exp2, avTestLum - avStudyLum)

exp2$absLumDiff = abs(exp2$avLumDiff)

# model
exp2_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp2, testChange==1), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp2_m1_loo=loo(exp2_m1)
exp2_m1_waic=waic(exp2_m1)

# include luminance as a predictor
exp2_m2 = update(exp2_m1, .~. + absLumDiff)

exp2_m2_loo=loo(exp2_m2)
exp2_m2_waic=waic(exp2_m2)

compare(exp2_m1_loo, exp2_m2_loo)
compare(exp2_m1_waic, exp2_m2_waic)

# predictive accuracy (WAIC) is improved by having luminance difference in the model
summary(exp2_m2) # bigger luminance diff = more change of responding change

with(subset(exp2, testChange ==1), hist(absLumDiff)) # however the range of luminance differences was pretty restricted

# luminance differences change only slightly with number of changes
aggregate(absLumDiff ~ Cond + SetSize, data = subset(exp2,testChange==1), FUN = mean)
# it is hard to see how this can explain changes to false-alarm rate

## Experiment 3
exp3 <- read.table(file = 'data/exp3/exp3_full')

exp3 = within(exp3, {
  Cond = as.factor(Cond)
  SetSize = as.factor(SetSize)
  ppt = as.factor(ppt)
})

exp3$avStudyLum = extractMeanRelLum(colorVec = exp3$study, N = exp3$SetSize)
exp3$avTestLum = extractMeanRelLum(colorVec = exp3$test, N = rep(1, length(exp3$test)))

exp3$avLumDiff = with(exp3, avTestLum - avStudyLum)

exp3$absLumDiff = abs(exp3$avLumDiff)

aggregate(exp3$respChange ~ cut(exp3$absLumDiff, 7)+exp3$testChange, FUN=mean)

# model
# change trials
exp3c_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp3, testChange==1), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp3c_m1_loo=loo(exp3c_m1)
exp3c_m1_waic=waic(exp3c_m1)

# include luminance as a predictor
exp3c_m2 = update(exp3c_m1, .~. + absLumDiff)

exp3c_m2_loo=loo(exp3c_m2)
exp3c_m2_waic=waic(exp3c_m2)

compare(exp3c_m1_loo, exp3c_m2_loo)
compare(exp3c_m1_waic, exp3c_m2_waic)

# including doesn't improve predictive accuracy for change trials
summary(exp3c_m2) # coefficient in the wrong direction

with(subset(exp3, testChange ==1), hist(absLumDiff))

# same trials
exp3s_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp3, testChange==0), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp3s_m1_loo=loo(exp3s_m1)
exp3s_m1_waic=waic(exp3s_m1)

# include luminance as a predictor
exp3s_m2 = update(exp3s_m1, .~. + absLumDiff)

exp3s_m2_loo=loo(exp3s_m2)
exp3s_m2_waic=waic(exp3s_m2)

compare(exp3s_m1_loo, exp3s_m2_loo)
compare(exp3s_m1_waic, exp3s_m2_waic)

# smaller luminance differences reduce the chance of responding change for same trials
summary(exp3s_m2)

with(subset(exp3, testChange == 0), hist(absLumDiff))

## Experiment 4
exp4 <- read.table(file = 'data/exp4/exp4_full')

exp4 = within(exp4, {
  Cond = as.factor(Cond)
  SetSize = as.factor(SetSize)
  ppt = as.factor(ppt)
})

exp4$avStudyLum = extractMeanRelLum(colorVec = exp4$study, N = exp4$SetSize)
exp4$avTestLum = extractMeanRelLum(colorVec = exp4$test, N = rep(1, length(exp4$test)))

exp4$avLumDiff = with(exp4, avTestLum - avStudyLum)

exp4$absLumDiff = abs(exp4$avLumDiff)

# model
# change trials
exp4c_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp4, testChange==1), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp4c_m1_loo=loo(exp4c_m1)
exp4c_m1_waic=waic(exp4c_m1)

# include luminance as a predictor
exp4c_m2 = update(exp4c_m1, .~. + absLumDiff)

exp4c_m2_loo=loo(exp4c_m2)
exp4c_m2_waic=waic(exp4c_m2)

compare(exp4c_m1_loo, exp4c_m2_loo)
compare(exp4c_m1_waic, exp4c_m2_waic)

# no effect of luminance for change trials
summary(exp4c_m2) # small coefficient (given range of differences)

with(subset(exp4, testChange == 1), hist(absLumDiff))

# same trials
exp4s_m1 = stan_glmer(respChange ~ 1 + Cond + SetSize + (1 | ppt), data = subset(exp4, testChange==0), family = binomial(link = logit), prior_intercept = cauchy(0, 2.5), prior = cauchy(0, 2.5))

exp4s_m1_loo=loo(exp4s_m1)
exp4s_m1_waic=waic(exp4s_m1)

# include luminance as a predictor
exp4s_m2 = update(exp4s_m1, .~. + absLumDiff)

exp4s_m2_loo=loo(exp4s_m2)
exp4s_m2_waic=waic(exp4s_m2)

compare(exp4s_m1_loo, exp4s_m2_loo)
compare(exp4s_m1_waic, exp4s_m2_waic)

# again, less probability of responding change to same trials with a small luminance difference between study and test
summary(exp4s_m2)

with(subset(exp4, testChange == 0), hist(absLumDiff))

save.image('lumAnalysis.RData')
