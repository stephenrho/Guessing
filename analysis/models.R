# JAGS models for analysis of experiments in 
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

setwd('analysis/')
dir.create('models/')
setwd('models/')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DATA MODEL --------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dataModel <- "
  model{
  # estimate hit and false alarm rates for experiments 1, 2, 3, and 4
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      logit(P[i]) <- Rsubj[testLevel[i], condLevel[i], ssLevel[i], id[i]] # freely estimate rates for each condition (with shrinkage)
    }

    # participant level parameters
    for (t in 1:2){ # same, change
      for (l in 1:N_c){ # condition (base rate manipulation or number of changes)
        for (ss in 1:N_ss){ # set sizes
          for (s in 1:S){
            Rsubj[t,l,ss,s] ~ dnorm(R[t,l,ss], R_Tau)
          }
        }
      }
    }

    # hierarchical priors
    # prior for means
    for (t in 1:2){
      for (l in 1:N_c){
        for (ss in 1:N_ss){
          R[t,l,ss] ~ dnorm(0, 1/10^2)
        }
      }
    }
    # standard deviations
    R_Tau <- 1/pow(R_SD, 2)
    R_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(text = dataModel, con = 'dataModel.txt')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# WHOLE DISPLAY MODELS ----------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

informedWD <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*g[i]) + (1 - a[i])*u[i], # p(hit)
                     a[i]*g[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # informed guessing
      g[i] <- ((1 - d[i])*u[i])/ ((1 - d[i])*u[i] + (1 - u[i]))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(informedWD, "informedWD.txt")

# -----------------------------------------------------------------------

informedWD_vk <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*g[i]) + (1 - a[i])*u[i], # p(hit)
                     a[i]*g[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # informed guessing
      g[i] <- ((1 - d[i])*u[i])/ ((1 - d[i])*u[i] + (1 - u[i]))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[ssLevel[i], id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    
    # participant level parameters
    for (ss in 1:N_ss){
      for (s in 1:S){
        Ksubj[ss, s] ~ dnorm(K[ss], K_Tau)T(,8)
      }
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # hierarchical priors
    # grand means
    K[1] ~ dnorm(3, 1/K1_sd^2) # to stabilize estimate at SS2 allow SD to be specified

    for (ss in 2:N_ss){
      K[ss] ~ dnorm(3, 1/10^2)
    }
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(informedWD_vk, "informedWD_vk.txt")

# -----------------------------------------------------------------------

uninformedWD <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*u[i]) + (1 - a[i])*u[i], # p(hit)
                     ifelse(d[i] == 1, (1 - a[i])*u[i], u[i])) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedWD, "uninformedWD.txt")

# -----------------------------------------------------------------------

uninformedWD_vk <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*u[i]) + (1 - a[i])*u[i], # p(hit)
                     ifelse(d[i] == 1, (1 - a[i])*u[i], u[i])) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[ssLevel[i], id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    # participant level parameters
    for (ss in 1:N_ss){
      for (s in 1:S){
        Ksubj[ss, s] ~ dnorm(K[ss], K_Tau)T(,8)
      }
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # grand means
    K[1] ~ dnorm(3, 1/K1_sd^2) # to stabilize estimate at SS2 allow SD to be specified

    for (ss in 2:N_ss){
      K[ss] ~ dnorm(3, 1/10^2)
    }
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedWD_vk, "uninformedWD_vk.txt")

# -----------------------------------------------------------------------

uninformedWD_vu <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*u[i]) + (1 - a[i])*u[i], # p(hit)
                     ifelse(d[i] == 1, (1 - a[i])*u[i], u[i])) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[ssLevel[i], B_level[i], id[i]] # separate u sample for each set size x base rate x participant
    }
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (ss in 1:N_ss){
      for (l in 1:B_n){
        for (s in 1:S){
          Usubj[ss, l, s] ~ dnorm(U[ss, l], U_Tau) # levels share variance
        }
      }
    }

    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (ss in 1:N_ss){
      for (l in 1:B_n){
        U[ss, l] ~ dnorm(0, 1/10^2)
      }
    }

    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedWD_vu, "uninformedWD_vu.txt")

# -----------------------------------------------------------------------

mixtureWD <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1,
                     a[i]*(d[i] + (1 - d[i])*lambda[i]) + (1 - a[i])*u[i], # p(hit)
                     ifelse(d[i] == 1, (1 - a[i])*u[i], a[i]*lambda[i] + (1 - a[i])*u[i]))# p(false-alarm)
      
      # mixture of guessing types
      lambda[i] <- (P_og[i]*g_o[i] + (1 - P_og[i])*g_i[i])
      
      # optimal guessing
      g_o[i] <- ifelse(g_i[i] == 0.5, 0.5, ifelse(g_i[i] > 0.5, 1, 0))
      
      # informed guessing
      g_i[i] <- ((1 - d[i])*u[i])/ ((1 - d[i])*u[i] + (1 - u[i]))
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]]
      logit(P_og[i]) <- P_og_subj[id[i]]
    }
    # individual level
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # Hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
    
    # mixture parameter
    for (s in 1:S){
      P_og_subj[s] ~ dnorm(POG, POG_Tau)
    }
    
    # priors
    POG ~ dnorm(0, 1/10^2)

    POG_Tau <- 1/pow(POG_SD, 2)
    POG_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(mixtureWD, "mixtureWD.txt")

# -----------------------------------------------------------------------

logisticWD <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*(d[i] + (1 - d[i])*o[i]) + (1 - a[i])*u[i], # p(hit)
                     a[i]*o[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      ## log factorial method for implementing choose function
      d[i] <- ifelse(k[i] >= N[i] - C[i] + 1, 1, 1 - exp(logfact(N[i] - k[i]) - (logfact(C[i]) + logfact(N[i] - k[i] - C[i]))) / exp(logfact(N[i]) - (logfact(C[i]) + logfact(N[i] - C[i]))))
      
      # informed guessing
      g[i] <- ((1 - d[i])*u[i])/ ((1 - d[i])*u[i] + (1 - u[i]))
      
      # logistic rule
      o[i] <- 1/(1 + exp(-lambda[i]*log(g[i]/(1 - g[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
      lambda[i] <- exp(Lsubj[id[i]])
    }
    
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
      Asubj[s] ~ dnorm(A, A_Tau)
      Lsubj[s] ~ dnorm(L, L_Tau)
      for (l in 1:B_n){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }

    # hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    L ~ dnorm(0, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
    L_Tau <- 1/pow(L_SD, 2)
    L_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(logisticWD, "logisticWD.txt")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# SINGLE PROBE MODELS -----------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

informedSP <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*g[i] + (1 - a[i])*u[i], # p(hit)
                     a[i]*(1 - d[i])*g[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # informed guessing
      g[i] <- u[i]/(u[i] + (1 - d[i])*(1 - u[i]))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(informedSP, "informedSP.txt")

# -----------------------------------------------------------------------

informedSP_vk <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*g[i] + (1 - a[i])*u[i], # p(hit)
                     a[i]*(1 - d[i])*g[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # informed guessing
      g[i] <- u[i]/(u[i] + (1 - d[i])*(1 - u[i]))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[ssLevel[i], id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    
    # participant level parameters
    for (ss in 1:N_ss){
      for (s in 1:S){
        Ksubj[ss, s] ~ dnorm(K[ss], K_Tau)
      }
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # hierarchical priors
    # grand means
    K[1] ~ dnorm(3, 1/K1_sd^2) # to stabilize estimate at SS2 allow SD to be specified

    for (ss in 2:N_ss){
      K[ss] ~ dnorm(3, 1/10^2)
    }
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(informedSP_vk, "informedSP_vk.txt")

# -----------------------------------------------------------------------

uninformedSP <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1,
                     ifelse(d[i] == 1, 1 - (1 - a[i])*(1 - u[i]), u[i]), # p(hit)
                     a[i]*(1 - d[i])*u[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)

      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedSP, "uninformedSP.txt")

# -----------------------------------------------------------------------

uninformedSP_vk <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1,
                     ifelse(d[i] == 1, 1 - (1 - a[i])*(1 - u[i]), u[i]), # p(hit)
                     a[i]*(1 - d[i])*u[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[ssLevel[i], id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
    }
    # participant level parameters
    for (ss in 1:N_ss){
      for (s in 1:S){
        Ksubj[ss, s] ~ dnorm(K[ss], K_Tau)
      }
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # grand means
    K[1] ~ dnorm(3, 1/K1_sd^2) # to stabilize estimate at SS2 allow SD to be specified

    for (ss in 2:N_ss){
      K[ss] ~ dnorm(3, 1/10^2)
    }
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedSP_vk, "uninformedSP_vk.txt")

# -----------------------------------------------------------------------

uninformedSP_vu <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1,
                     ifelse(d[i] == 1, 1 - (1 - a[i])*(1 - u[i]), u[i]), # p(hit)
                     a[i]*(1 - d[i])*u[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[ssLevel[i], B_level[i], id[i]] # separate u sample for each set size x base rate x participant
    }
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
    }
    for (s in 1:S){
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (ss in 1:N_ss){
      for (l in 1:B_n){
        for (s in 1:S){
          Usubj[ss, l, s] ~ dnorm(U[ss, l], U_Tau) # levels share variance
        }
      }
    }

    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (ss in 1:N_ss){
      for (l in 1:B_n){
        U[ss, l] ~ dnorm(0, 1/10^2)
      }
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(uninformedSP_vu, "uninformedSP_vu.txt")

# -----------------------------------------------------------------------

mixtureSP <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1,
                     ifelse(d[i]==1, a[i] + (1-a[i])*u[i], a[i]*lambda[i] + (1-a[i])*u[i]), # p(hit)
                     a[i]*(1 - d[i])*lambda[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      # mixture of guessing types
      lambda[i] <- (P_og[i]*g_o[i] + (1 - P_og[i])*g_i[i])
      
      # optimal guessing
      g_o[i] <- ifelse(g_i[i] == 0.5, 0.5, ifelse(g_i[i] > 0.5, 1, 0))
      
      # informed guessing
      g_i[i] <- u[i]/(u[i] + (1 - d[i])*(1 - u[i]))
      
      d[i] <- min(k[i]/N[i], 1)
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]]
      logit(P_og[i]) <- P_og_subj[id[i]]
    }
    # individual level
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
      Asubj[s] ~ dnorm(A, A_Tau)
    }
    for (l in 1:B_n){
      for (s in 1:S){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # Hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
    
    # mixture parameter
    for (s in 1:S){
      P_og_subj[s] ~ dnorm(POG, POG_Tau)
    }
    
    # priors
    POG ~ dnorm(0, 1/10^2)

    POG_Tau <- 1/pow(POG_SD, 2)
    POG_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(mixtureSP, "mixtureSP.txt")

# -----------------------------------------------------------------------

logisticSP <- "
  model{
    for (i in 1:n){
      #y[i] ~ dbern(y.hat[i])
      y[i] ~ dbin(y.hat[i], N_trials[i])
      y.hat[i] <- max(0, min(1, P[i]))
      
      # probability of responding change
      P[i] <- ifelse(CHANGE[i] == 1, 
                     a[i]*o[i] + (1 - a[i])*u[i], # p(hit)
                     a[i]*(1 - d[i])*o[i] + (1 - a[i])*u[i]) # p(false-alarm)
      
      d[i] <- min(k[i]/N[i], 1)
      
      # informed guessing
      g[i] <- u[i]/(u[i] + (1 - d[i])*(1 - u[i]))
      
      # logistic rule
      o[i] <- 1/(1 + exp(-lambda[i]*log(g[i]/(1 - g[i]))))
      
      # model transformations of k, u, and a
      k[i] <- max(kappa[i], 0) # Mass-at-chance transformation
      kappa[i] <- Ksubj[id[i]]
      logit(a[i]) <- Asubj[id[i]]
      logit(u[i]) <- Usubj[B_level[i], id[i]] # separate u sample for each base rate x participant
      lambda[i] <- exp(Lsubj[id[i]])
    }
    
    # participant level parameters
    for (s in 1:S){
      Ksubj[s] ~ dnorm(K, K_Tau)
      Asubj[s] ~ dnorm(A, A_Tau)
      Lsubj[s] ~ dnorm(L, L_Tau)
      for (l in 1:B_n){
        Usubj[l, s] ~ dnorm(U[l], U_Tau) # levels share variance
      }
    }
    
    # hierarchical priors
    # grand means
    K ~ dnorm(3, 1/10^2)
    A ~ dnorm(3, 1/10^2)
    L ~ dnorm(0, 1/10^2)
    for (l in 1:B_n){
      U[l] ~ dnorm(0, 1/10^2)
    }
    
    # standard deviations
    K_Tau <- 1/pow(K_SD, 2)
    K_SD ~ dgamma(1.01005, 0.1005012) # mode = .1, SD = 10 (v. vauge)
    U_Tau <- 1/pow(U_SD, 2)
    U_SD ~ dgamma(1.01005, 0.1005012)
    A_Tau <- 1/pow(A_SD, 2)
    A_SD ~ dgamma(1.01005, 0.1005012)
    L_Tau <- 1/pow(L_SD, 2)
    L_SD ~ dgamma(1.01005, 0.1005012)
  }"

writeLines(logisticSP, "logisticSP.txt")
