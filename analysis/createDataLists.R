### RE-FORMAT DATA INTO LISTS FOR USE WITH JAGS
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

rm(list=ls())

library(plyr)

# READ DATA ----------------------------------------------------
setwd('analysis/')

## for each experiment:
# read in raw data
# collapse trial level data counting number of trials and number of 'change' responses

#### EXPERIMENT 1 ------
exp1 <- read.table(file = 'data/exp1/exp1_full')
exp1$Cond <- as.factor(exp1$Cond)

exp1_coll = ddply(exp1, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
      N_respChange = sum(respChange),
      N_trials = length(respChange))

DATALIST_EXP1 <- list(
  y = exp1_coll$N_respChange,
  N_trials = exp1_coll$N_trials,
  CHANGE = exp1_coll$testChange,
  C = rep(1, nrow(exp1_coll)),
  B_level = as.numeric(exp1_coll$Cond),
  B_n = length(unique(exp1_coll$Cond)),
  N = exp1_coll$SetSize,
  S = length(unique(exp1_coll$ppt)),
  id = exp1_coll$ppt,
  n = nrow(exp1_coll),
  # for data model
  condLevel = as.numeric(exp1_coll$Cond),
  N_c = length(unique(exp1_coll$Cond)),
  ssLevel = as.numeric(as.factor(exp1_coll$SetSize)),
  N_ss = length(unique(exp1_coll$SetSize)),
  testLevel = as.numeric(as.factor(exp1_coll$testChange))
)

#### EXPERIMENT 2 ------
exp2 <- read.table(file = 'data/exp2/exp2_full')

exp2_coll = ddply(exp2, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  N_respChange = sum(respChange),
                  N_trials = length(respChange))

DATALIST_EXP2 <- list(
  y = exp2_coll$N_respChange,
  N_trials = exp2_coll$N_trials,
  CHANGE = exp2_coll$testChange,
  C = exp2_coll$Cond,
  B_level = rep(1, nrow(exp2_coll)),
  B_n = 1,
  N = exp2_coll$SetSize,
  S = length(unique(exp2_coll$ppt)),
  id = exp2_coll$ppt,
  n = nrow(exp2_coll),
  # for data model
  condLevel = exp2_coll$Cond,
  N_c = length(unique(exp2_coll$Cond)),
  ssLevel = as.numeric(as.factor(exp2_coll$SetSize)),
  N_ss = length(unique(exp2_coll$SetSize)),
  testLevel = as.numeric(as.factor(exp2_coll$testChange))
)

#### EXPERIMENT 3 ------
exp3 <- read.table(file = 'data/exp3/exp3_full')
exp3$Cond <- as.factor(exp3$Cond)

exp3_coll = ddply(exp3, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  N_respChange = sum(respChange),
                  N_trials = length(respChange))

DATALIST_EXP3 <- list(
  y = exp3_coll$N_respChange,
  N_trials = exp3_coll$N_trials,
  CHANGE = exp3_coll$testChange,
  B_level = as.numeric(exp3_coll$Cond),
  B_n = length(unique(exp3_coll$Cond)),
  N = exp3_coll$SetSize,
  S = length(unique(exp3_coll$ppt)),
  id = exp3_coll$ppt,
  n = nrow(exp3_coll),
  # for data model
  condLevel = as.numeric(exp3_coll$Cond),
  N_c = length(unique(exp3_coll$Cond)),
  ssLevel = as.numeric(as.factor(exp3_coll$SetSize)),
  N_ss = length(unique(exp3_coll$SetSize)),
  testLevel = as.numeric(as.factor(exp3_coll$testChange))
)

#### EXPERIMENT 4 ------
exp4 <- read.table(file = 'data/exp4/exp4_full')
exp4$Cond <- as.factor(exp4$Cond)

exp4_coll = ddply(exp4, c('ppt', 'SetSize', 'Cond', 'testChange'), summarize,
                  N_respChange = sum(respChange),
                  N_trials = length(respChange))

DATALIST_EXP4 <- list(
  y = exp4_coll$N_respChange,
  N_trials = exp4_coll$N_trials,
  CHANGE = exp4_coll$testChange,
  B_level = as.numeric(exp4_coll$Cond),
  B_n = length(unique(exp4_coll$Cond)),
  N = exp4_coll$SetSize,
  S = length(unique(exp4_coll$ppt)),
  id = exp4_coll$ppt,
  n = nrow(exp4_coll),
  # for data model
  condLevel = as.numeric(exp4_coll$Cond),
  N_c = length(unique(exp4_coll$Cond)),
  ssLevel = as.numeric(as.factor(exp4_coll$SetSize)),
  N_ss = length(unique(exp4_coll$SetSize)),
  testLevel = as.numeric(as.factor(exp4_coll$testChange))
)

#### SAVE DATA LISTS
save(list = ls(pattern = 'DATALIST'), file = 'guessingDataLists.RData')

