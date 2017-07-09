# functions for creating vectors of likelihoods at each step of the chain for a given observation
# functions used in calculation of WAIC
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

lik_dataM <- function(dataList, obs, samples, predFunc=NULL, logLik = F){
  
  # [testLevel, condLevel, ssLevel, id]
  
  # get observation
  y = dataList$y[obs]
  N_trials = dataList$N_trials[obs]
  id = dataList$id[obs]
  testLevel = dataList$testLevel[obs]
  condLevel = dataList$condLevel[obs]
  ssLevel = dataList$ssLevel[obs]
  id = dataList$id[obs]
  
  pcol = paste0('Rsubj[', testLevel, ',', condLevel,',', ssLevel,',', id, ']' )
  
  lik = rep(0, nrow(samples)) # initialize vector of likelihoods
  for (i in 1:nrow(samples)){
    step = samples[i,]
    
    p = plogis(step[[pcol]])
    lik[i] = dbinom(y, N_trials, p, log = logLik)
  }
  return(lik)
}

lik_guessing <- function(dataList, obs, samples, predFunc, logLik = F){
  # works out likelihood of an observation at each step in chain given hierarchical parameters
  # for standard guessing models
  
  # get observation
  y = dataList$y[obs]
  N_trials = dataList$N_trials[obs]
  id = dataList$id[obs]
  CHANGE = dataList$CHANGE[obs]
  B_level = dataList$B_level[obs]
  N = dataList$N[obs]
  
  if ("C" %in% names(dataList)){ # allows for SP data
    C = dataList$C[obs]
  } else{
    C = NULL
  }
  
  if ('K[1]' %in% colnames(samples)){ # allows for the variable k version
    kcol = paste0('Ksubj[', B_level, ',', id, ']')
  } else{
    kcol = paste0('Ksubj[', id, ']')
  }
  
  lik = rep(0, nrow(samples)) # initialize vector of likelihoods
  for (i in 1:nrow(samples)){
    step = samples[i,]
    
    # sample from hyperparameters
    a = plogis(step[[paste0('Asubj[', id, ']')]])
    k = max(0 , step[[kcol]])
    u = plogis(step[[paste0('Usubj[', B_level, ',', id, ']')]])
    
    # likelihood
    if (is.null(C)){
      p = predFunc(k = k, a = a, u = u, n = N)
    } else{
      p = predFunc(k = k, a = a, u = u, n = N, nc = C)
    }
    
    p = ifelse(CHANGE == 1, p[1], p[3])
    
    lik[i] = dbinom(y, N_trials, p, log = logLik)
  }
  return(lik)
}

lik_mixture <- function(dataList, obs, samples, predFunc, logLik = F){
  # works out likelihood of an observation at each step in chain given hierarchical parameters
  # for mixture model (additional guessing model)
  
  # get observation
  y = dataList$y[obs]
  N_trials = dataList$N_trials[obs]
  id = dataList$id[obs]
  CHANGE = dataList$CHANGE[obs]
  B_level = dataList$B_level[obs]
  N = dataList$N[obs]
  
  if ("C" %in% names(dataList)){ # allows for SP data
    C = dataList$C[obs]
  } else{
    C = NULL
  }
  
  lik = rep(0, nrow(samples)) # initialize vector of likelihoods
  for (i in 1:nrow(samples)){
    step = samples[i,]
    
    # sample from hyperparameters
    a = plogis(step[[paste0('Asubj[', id, ']')]])
    k = max(0 , step[[paste0('Ksubj[', id, ']')]])
    u = plogis(step[[paste0('Usubj[', B_level, ',', id, ']')]])
    pog = plogis(step[[paste0('P_og_subj[', id, ']')]])
    
    # likelihood
    if (is.null(C)){
      p = predFunc(k = k, a = a, u = u, n = N, P_og = pog)
    } else{
      p = predFunc(k = k, a = a, u = u, n = N, nc = C, P_og = pog)
    }
    
    p = ifelse(CHANGE == 1, p[1], p[3])
    
    lik[i] = dbinom(y, N_trials, p, log = logLik)
  }
  return(lik)
}

lik_logistic <- function(dataList, obs, samples, predFunc, logLik = F){
  # works out likelihood of an observation at each step in chain given hierarchical parameters
  # for logistic model (additional guessing model)
  
  # get observation
  y = dataList$y[obs]
  N_trials = dataList$N_trials[obs]
  id = dataList$id[obs]
  CHANGE = dataList$CHANGE[obs]
  B_level = dataList$B_level[obs]
  N = dataList$N[obs]
  
  if ("C" %in% names(dataList)){ # allows for SP data
    C = dataList$C[obs]
  } else{
    C = NULL
  }
  
  lik = rep(0, nrow(samples)) # initialize vector of likelihoods
  for (i in 1:nrow(samples)){
    step = samples[i,]
    
    # sample from hyperparameters
    a = plogis(step[[paste0('Asubj[', id, ']')]])
    k = max(0 , step[[paste0('Ksubj[', id, ']')]])
    u = plogis(step[[paste0('Usubj[', B_level, ',', id, ']')]])
    l = exp(step[[paste0('Lsubj[', id, ']')]])
    
    # likelihood
    if (is.null(C)){
      p = predFunc(k = k, a = a, u = u, n = N, l = l)
    } else{
      p = predFunc(k = k, a = a, u = u, n = N, nc = C, l = l)
    }
    
    p = ifelse(CHANGE == 1, p[1], p[3])
    
    lik[i] = dbinom(y, N_trials, p, log = logLik)
  }
  return(lik)
}

# a master function used to calculate WAIC
getWAIC <- function(dataList, samples, likFun, predFunc){
  # see http://kylehardman.com/BlogPosts/View/7
  # each dataList will have an element y
  
  nObvs = length(dataList$y)
  
  # totals
  lppd = 0
  p_1 = 0
  p_2 = 0
  
  for (y_i in 1:nObvs){
    lik = likFun(dataList = dataList, obs = y_i, samples = samples, predFunc=predFunc)
    
    logMeanLik = log(mean(lik)) # used twice
    
    lppd = lppd + logMeanLik
    
    # two estimates of effective number of parameters
    p_1 = p_1 + 2*(logMeanLik - mean(log(lik)))
    
    p_2 = p_2 + var(log(lik))
  }
  # calculate the two waic versions (on deviance scale)
  waic_1 = -2*(lppd - p_1)
  waic_2 = -2*(lppd - p_2)
  
  return(list(lppd=lppd, p_1=p_1, p_2=p_2, waic_1=waic_1, waic_2=waic_2))
}