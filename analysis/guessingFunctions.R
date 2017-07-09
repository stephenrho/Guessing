#### FUNCTIONS USED DURING ANALYSIS AND WRITE UP
# Rhodes et al. - Informed Guessing in Change Detection
# Author: Stephen Rhodes (rhodessp at missouri.edu)
# License: GNU GPL v3.0

#### EXPERIMENT 1 ------

plotDataM_1 <- function(add = T, errWidth = .01){
  if (!add){
    plot(1000, xlim=c(0,1), ylim=c(0,1), xlab='', ylab='', main='', axes=F)
    box()
  }
  N = c(2,5,8)
  P = c(0.3, 0.5, 0.7)
  
  for (n in N){
    for (p in P){
      # data to plot 
      toplot <- subset(fh.estimates.postE1, N == n & P == p)
      # ERROR BARS
      # f
      with(toplot, segments(x0 = f.lower, x1 = f.upper,
                            y0 = h.median, y1 = h.median))
      with(toplot, segments(x0 = f.lower, x1 = f.lower,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      with(toplot, segments(x0 = f.upper, x1 = f.upper,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      # h
      with(toplot, segments(y0 = h.lower, y1 = h.upper,
                            x0 = f.median, x1 = f.median))
      with(toplot, segments(y0 = h.lower, y1 = h.lower,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      with(toplot, segments(y0 = h.upper, y1 = h.upper,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      # (f,h) point
      with(toplot, points(f.median, h.median, pch = 16, col = 'black', cex = .7))
    }
  } 
}

# model functions
uninformed_1 <- function(k,a,u,n){
  
  d = min(k/n, 1)
  
  p = 1:4
  p[1] = a*(d + (1 - d)*u) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = ifelse(d == 1, (1 - a)*u, u)
  p[4] = 1 - p[3]
  return(p)
}

informed_1 <- function(k,a,u,n){
  
  d = min(k/n, 1)
  
  g = ((1 - d)*u)/((1 - d)*u + (1 - u))
  
  p = 1:4
  p[1] = a*(d + (1 - d)*g) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*g + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

mixture_1 <- function(k,a,u,n,P_og){
  
  d = min(k/n, 1)
  
  g_i = ((1 - d)*u)/((1 - d)*u + (1 - u))
  g_o = ifelse(g_i == 0.5, 0.5, ifelse(g_i > 0.5, 1, 0))
  
  lambda <- (P_og*g_o + (1 - P_og)*g_i)
  
  p = 1:4
  p[1] = a*(d + (1 - d)*lambda) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = ifelse(d == 1, (1 - a)*u, a*lambda + (1 - a)*u)
  p[4] = 1 - p[3]
  return(p)
}

logistic_1 <- function(k,a,u,n,l){
  
  d = min(k/n, 1)
  
  g <- ((1 - d)*u)/ ((1 - d)*u + (1 - u))
  
  o <- 1/(1 + exp(-l*log(g/(1 - g))))
  
  p = 1:4
  p[1] = a*(d + (1 - d)*o) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*o + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

# use the above functions to generate (f,h) points for each posterior sample
genPredsInformedExp1 <- function(k, a, u, long=F){
  Ns = c(2, 5, 8); Ps = c(0.3, 0.5, 0.7)
  tmp <- c()
  for (p in 1:length(Ps)){
    for (n in 1:length(Ns)){
      fh = informed_1(k = k, a = a, u = u[p], n = Ns[n])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], P = Ps[p], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'P'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

genPredsUninformedExp1 <- function(k, a, u, long=F){
  Ns = c(2, 5, 8); Ps = c(0.3, 0.5, 0.7)
  tmp <- c()
  for (p in 1:length(Ps)){
    for (n in 1:length(Ns)){
      fh = uninformed_1(k = k, a = a, u = u[p], n = Ns[n])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], P = Ps[p], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'P'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

#### EXPERIMENT 2 -----

plotDataM_2 <- function(add = T, errWidth = .1){
  if (!add){
    plot(1000, xlim=c(0.5,9.5), ylim=c(0,1), xlab='', ylab='', main='', axes=F)
    box()
  }
  N = c(5,8)
  
  for (n in N){
    toplot <- subset(fh.estimates.postE2, N == n)
    
    # determine x locations
    if (n == 5){
      x = 1:4
    } else{
      x = 6:9
    }
    
    # f
    with(toplot, segments(x0 = x, x1 = x, y0 = f.lower, y1 = f.upper))
    with(toplot, segments(x0 = x - errWidth, x1 = x + errWidth, y0 = f.lower, y1 = f.lower))
    with(toplot, segments(x0 = x - errWidth, x1 = x + errWidth, y0 = f.upper, y1 = f.upper))
    
    with(toplot, points(x, f.median, pch = 16, col = 'black', cex = .7, type='b'))
    
    # h
    with(toplot, segments(x0 = x, x1 = x, y0 = h.lower, y1 = h.upper))
    with(toplot, segments(x0 = x - errWidth, x1 = x + errWidth, y0 = h.lower, y1 = h.lower))
    with(toplot, segments(x0 = x - errWidth, x1 = x + errWidth, y0 = h.upper, y1 = h.upper))
    
    with(toplot, points(x, h.median, pch = 16, col = 'black', cex = .7, type='b'))
  }
}

# Model functions
uninformed_2 <- function(k,a,u,n,nc){
  
  d = ifelse(k >= n - nc + 1, 1, 1 - choose(n - k, nc)/choose(n, nc))
  
  p = 1:4
  p[1] = a*(d + (1 - d)*u) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = ifelse(d == 1, (1 - a)*u, u)
  p[4] = 1 - p[3]
  return(p)
}

informed_2 <- function(k,a,u,n,nc){
  
  d = ifelse(k >= n - nc + 1, 1, 1 - choose(n - k, nc)/choose(n, nc))
  
  g = ((1 - d)*u)/((1 - d)*u + (1 - u))
  
  p = 1:4
  p[1] = a*(d + (1 - d)*g) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*g + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

mixture_2 <- function(k,a,u,n,nc,P_og){
  
  d = ifelse(k >= n - nc + 1, 1, 1 - choose(n - k, nc)/choose(n, nc))
  
  g_i = ((1 - d)*u)/((1 - d)*u + (1 - u))
  g_o = ifelse(g_i == 0.5, 0.5, ifelse(g_i > 0.5, 1, 0))
  
  lambda <- (P_og*g_o + (1 - P_og)*g_i)
  
  p = 1:4
  p[1] = a*(d + (1 - d)*lambda) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = ifelse(d == 1, (1 - a)*u, a*lambda + (1 - a)*u)
  p[4] = 1 - p[3]
  return(p)
}

logistic_2 <- function(k,a,u,n,nc,l){
  
  d = ifelse(k >= n - nc + 1, 1, 1 - choose(n - k, nc)/choose(n, nc))
  
  g <- ((1 - d)*u)/ ((1 - d)*u + (1 - u))
  
  o <- 1/(1 + exp(-l*log(g/(1 - g))))
  
  p = 1:4
  p[1] = a*(d + (1 - d)*o) + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*o + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

# Generate predictions for specific versions of the models

genPredsInformedExp2 <- function(k, a, u, long=F){
  Ns = c(5, 8); Cs = 1:4
  tmp <- c()
  for (n in 1:length(Ns)){
    for (nc in 1:length(Cs)){
      fh = informed_2(k = k, a = a, u = u, n = Ns[n], nc = Cs[nc])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], C = Cs[nc], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(5,8), each = 2*4), rep(1:4, each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'C'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

genPredsUninformedExp2 <- function(k, a, u, long=F){
  Ns = c(5, 8); Cs = 1:4
  tmp <- c()
  for (n in 1:length(Ns)){
    for (nc in 1:length(Cs)){
      fh = uninformed_2(k = k, a = a, u = u, n = Ns[n], nc = Cs[nc])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], C = Cs[nc], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(5,8), each = 2*4), rep(1:4, each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'C'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

#### EXPERIMENT 3 -----

plotDataM_3 <- function(add = T, errWidth = .01){
  if (!add){
    plot(1000, xlim=c(0,1), ylim=c(0,1), xlab='', ylab='', main='', axes=F)
    box()
  }
  N = c(2,5,8)
  P = c(0.3, 0.5, 0.7)
  
  for (n in N){
    for (p in P){
      # data to plot 
      toplot <- subset(fh.estimates.postE3, N == n & P == p)
      # ERROR BARS
      # f
      with(toplot, segments(x0 = f.lower, x1 = f.upper,
                            y0 = h.median, y1 = h.median))
      with(toplot, segments(x0 = f.lower, x1 = f.lower,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      with(toplot, segments(x0 = f.upper, x1 = f.upper,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      # h
      with(toplot, segments(y0 = h.lower, y1 = h.upper,
                            x0 = f.median, x1 = f.median))
      with(toplot, segments(y0 = h.lower, y1 = h.lower,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      with(toplot, segments(y0 = h.upper, y1 = h.upper,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      # (f,h) point
      with(toplot, points(f.median, h.median, pch = 16, col = 'black', cex = .7))
    }
  } 
}

# model functions
uninformed_3 <- function(k,a,u,n){
  
  d = min(k/n, 1)
  
  p = 1:4
  p[1] = ifelse(d == 1, 1 - (1 - a)*(1 - u), u)
  p[2] = 1 - p[1]
  p[3] = a*(1 - d)*u + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

informed_3 <- function(k,a,u,n){
  
  d = min(k/n, 1)
  
  g = u/(u + (1 - d)*(1 - u))
  
  p = 1:4
  p[1] = a*g + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*(1 - d)*g + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

mixture_3 <- function(k,a,u,n,P_og){
  
  d = min(k/n, 1)
  
  g_i =  u/(u + (1 - d)*(1 - u))
  g_o = ifelse(g_i == 0.5, 0.5, ifelse(g_i > 0.5, 1, 0))
  
  lambda <- (P_og*g_o + (1 - P_og)*g_i)
  
  p = 1:4
  p[1] = ifelse(d==1, a + (1-a)*u, a*lambda + (1-a)*u)
  p[2] = 1 - p[1]
  p[3] = a*(1 - d)*lambda + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

logistic_3 <- function(k,a,u,n,l){
  
  d = min(k/n, 1)
  
  g <- u/(u + (1 - d)*(1 - u))
  
  o <- 1/(1 + exp(-l*log(g/(1 - g))))
  
  p = 1:4
  p[1] = a*o + (1 - a)*u
  p[2] = 1 - p[1]
  p[3] = a*(1 - d)*o + (1 - a)*u
  p[4] = 1 - p[3]
  return(p)
}

# use the above generic models to generate predictions with specific parameters
genPredsInformedExp3 <- function(k, a, u, long=F){
  Ns = c(2, 5, 8); Ps = c(0.3, 0.5, 0.7)
  tmp <- c()
  for (p in 1:length(Ps)){
    for (n in 1:length(Ns)){
      fh = informed_3(k = k, a = a, u = u[p], n = Ns[n])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], P = Ps[p], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'P'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

genPredsUninformedExp3 <- function(k, a, u, long=F){
  Ns = c(2, 5, 8); Ps = c(0.3, 0.5, 0.7)
  tmp <- c()
  for (p in 1:length(Ps)){
    for (n in 1:length(Ns)){
      fh = uninformed_3(k = k, a = a, u = u[p], n = Ns[n])[c(3,1)]
      tmp <- rbind(tmp, c(N = Ns[n], P = Ps[p], f = fh[1], h = fh[2]))
    }
  }
  if (long){
    tmp2 <- cbind(rep(c(2,5,8), each = 2*3), rep(c(.3,.5,.7), each = 2), rep(c(0,1)), rep(0))
    for (r in 1:nrow(tmp2)){
      tmp2[r,4] <- tmp[tmp[,'N'] == tmp2[r,1] & tmp[,'P'] == tmp2[r,2], 3+tmp2[r,3]]
    }
    return(tmp2)
  } else{
    return(tmp)
  }
}

#### EXPERIMENT 4 -----

plotDataM_4 <- function(add = T, errWidth = .01){
  if (!add){
    plot(1000, xlim=c(0,1), ylim=c(0,1), xlab='', ylab='', main='', axes=F)
    box()
  }
  N = c(2,5,8)
  P = c(0.3, 0.5, 0.7)
  
  for (n in N){
    for (p in P){
      # data to plot 
      toplot <- subset(fh.estimates.postE4, N == n & P == p)
      # ERROR BARS
      # f
      with(toplot, segments(x0 = f.lower, x1 = f.upper,
                            y0 = h.median, y1 = h.median))
      with(toplot, segments(x0 = f.lower, x1 = f.lower,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      with(toplot, segments(x0 = f.upper, x1 = f.upper,
                            y0 = h.median - errWidth, y1 = h.median + errWidth))
      # h
      with(toplot, segments(y0 = h.lower, y1 = h.upper,
                            x0 = f.median, x1 = f.median))
      with(toplot, segments(y0 = h.lower, y1 = h.lower,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      with(toplot, segments(y0 = h.upper, y1 = h.upper,
                            x0 = f.median - errWidth, x1 = f.median + errWidth))
      # (f,h) point
      with(toplot, points(f.median, h.median, pch = 16, col = 'black', cex = .7))
    }
  } 
}

# other functions same as 3

#### OTHER USEFUL FUNCTIONS -----
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # from Kruschke (2015)
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

logistic <- function(x){
  1/(1 + exp(-x))
}

faintCol <- function(colLabs, alpha = 10){
  colVec = c()
  for (i in colLabs){
    RGB = col2rgb(i)
    newcol = rgb(red = RGB[1], green = RGB[2], blue = RGB[3], alpha = alpha, maxColorValue = 255)
    colVec = c(colVec, newcol)
  }
  return(colVec)
}

# useful abbrev function
printD <- function(num, decP){
  sprintf(fmt = paste0('%.', decP, 'f'), num)
}

getDIC <- function(samples, pD = F){
  dic = get(deparse(substitute(samples)))$BUGSoutput$DIC
  pd = get(deparse(substitute(samples)))$BUGSoutput$pD
  if (!pD){
    return(dic)
  } else{
    return(c(dic, pd))
  }
}

k_initFunc1 = function(){
  list('K' = runif(1, 0, 1))
       
}

k_initFunc2 = function(){
  initList = list()
  for (i in 1:2){
    initList[paste0('K[', i, ']')] = runif(1, 0, 1)
  }
  return(initList)
}

k_initFunc3 = function(){
  initList = list()
  for (i in 1:3){
    initList[paste0('K[', i, ']')] = runif(1, 0, 1)
  }
  return(initList)
}

