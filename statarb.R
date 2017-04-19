library(RcppFaddeeva)
library(sde)
library(doParallel)
library(foreach)

# setwd('C:/Users/Stavros/Google Drive/Trading/R/cointegration/')
# setwd('E:/R/statarb/R')
# source('util.R')
# source('PCA.R')

#' Fit OU model via Max Likelihood
#' @param S numeric input series
#' @param N length of input series
#' @param dt time step
#' @author Stav Sioutis
#' @rdname OUML
#' @export
OUML = function(S, N = length(S), dt = 1) 
{
  Sx = sum(S[1:(N-1)])
  Sy = sum(S[2:N])
  Sxx = sum(S[1:(N-1)]^2)
  Sxy = sum(S[1:(N-1)] * S[2:N])
  Syy = sum(S[2:N]^2)
  
  n = N-1
  mu = (Sy*Sxx - Sx*Sxy)/(n*(Sxx-Sxy) - (Sx^2 - Sx*Sy))
  lambda = -(1/dt)*log( (Sxy - mu*Sx - mu*Sy + n*mu^2) / (Sxx - 2*mu*Sx + n*mu^2));
  alpha = exp(-lambda*dt)
  alpha2 = exp(-2*lambda*dt)
  sigmahat = (1/n) * (Syy - 2*alpha*Sxy + alpha2*Sxx - 2*mu*(1-alpha)*(Sy - alpha*Sx)
                      +n * mu^2 * (1-alpha)^2)
  sigma = sqrt(sigmahat*2*lambda / (1-alpha2))
  
  return(list("mu" = mu, "sigma" = sigma, "lambda" = lambda))
}

#' Forecast OU model
#' @param S numeric input series
#' @param N length of input series
#' @param dt time step
#' @param color passthrough to \code{lines} plotting function
#' @param dWt Brownian motion
#' @author Stav Sioutis
#' @rdname OUforecast_ML
#' @export
OUforecast_ML = function(S, N = length(S), dt = 1, color = 'blue', dWt = sqrt(dt) * rnorm(N))
{
  vec = unlist(OUML(S, N, dt))
  
  mu = vec[1]
  sigma = vec[2]
  lambda = vec[3]
  
  cat(sigma, lambda, mu)
  
  # Simulate OU process
  S1 = rep(0, N)
  S1[1] = tail(S,1)
  
  for (t in 2:N) 
    S1[t] = S1[t-1] + lambda*(mu-S1[t-1])*dt + sigma*dWt[t]
  
  
  lines(S1, x=(1:N)+length(S),col=color)
}

#' helper function for calculating ET and VT
#' @param z numeric series
#' @param n Summation upper limit
#' @author Stav Sioutis
#' @rdname w1
w1 = function(z, n=50)
{
  k = 1:n
  w1.z = (0.5 * sum( gamma( k/2 )*( sqrt(2)*z )^k/factorial(k) ) )^2 - 
    (0.5 * sum( (-1)^k * gamma( k/2 )*( sqrt(2)*z )^k/factorial(k) ) )^2 
  return(w1.z)
}

#' helper function for calculating ET and VT
#' @param z numeric series
#' @param n Summation upper limit
#' @author Stav Sioutis
#' @rdname w2
w2 = function(z, n=50)
{
  k = 2*(1:n)-1
  w2.z = sum( gamma( k/2 ) * digamma(k/2)*( sqrt(2)*z)^k/factorial(k) ) 
  return(w2.z)
}

#' Returns the expected value of an entry point a given the parameters of an OU model
#' @param x entry point a. (Named x to satisfy optimization call)
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @param c transaction cost
#' @author Stav Sioutis
#' @rdname opt.a
opt.a = function(x, alpha, eta, c)
{
  a = x
  return(Re( exp(alpha*a^2 / eta^2)*(2*a + c) - eta * sqrt(pi/alpha) * erfi(a*sqrt(alpha)/eta)) )
}

#' Find trade entry point a which maximizes the expected value of the trade. Bertram (2011)
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @param c transaction cost
#' @author Stav Sioutis
#' @rdname opt.a
#' @return Optimal entry point a. Assumes optimal entry point m is symmetric about the mean of the OU model
#' @export
Max.Expected = function(alpha, c, eta)
{
  return(uniroot(opt.a, upper = 0, lower = -1, alpha = alpha, c = txn.cost, eta = eta)$root)
}

#' Expected value of trade. Bertram (2011)
#' @param a Trade entry point
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @author Stav Sioutis
#' @rdname ET
#' @export
ET = function(a, alpha, eta)
{
  m = -a
  return(( pi/alpha*(erfi( m*sqrt(alpha)/eta )-erfi( a*sqrt(alpha)/eta )) ))
}

#' Variance of a trade. Bertram (2011)
#' @param a Trade entry point
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @author Stav Sioutis
#' @rdname VT
#' @export
VT = function(a, alpha, eta)
{
  m = -a
  w.m = m*sqrt(2*alpha)/eta
  w.a = a*sqrt(2*alpha)/eta
  return( (w1(w.m) - w1(w.a) - w2(w.m) + w2(w.a))/alpha^2 )
}

#' Returns the expected Sharpe ratio of a trade entry point \code{a} given the parameters of an OU model
#' @param x entry point a. (Named x to satisfy optimization call)
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @param c transaction cost
#' @param rf risk-free rate
#' @author Stav Sioutis
#' @rdname opt.S
opt.S = function(x, alpha, eta, c = 0.01, rf = 0.001)
{
  a = x
  m = -a
  s.num = alpha * pi * (erfi(m*sqrt(alpha)/eta) - erfi(a*sqrt(alpha)/eta) )
  s.den = w1(m*sqrt(2*alpha)/eta) - w1(a*sqrt(2*alpha)/eta) - w2(m*sqrt(2*alpha)/eta) + w2(a*sqrt(2*alpha)/eta)
  S = (m-a-c-rf)/sqrt((m-a-c)^2) * sqrt(s.num / s.den)
  return(Re(S))
}

#' Find trade entry point that maximizes Sharpe ratio
#' @param alpha OU model parameter
#' @param eta OU model parameter
#' @param c transaction cost
#' @author Stav Sioutis
#' @rdname Max.Sharpe
#' @return Optimal entry point a. Assumes optimal entry point m is symmetric about the mean of the OU model
#' @export
Max.Sharpe = function(alpha, eta, c = 0.01, rf = 0.001)
{
  return(optimize(opt.S, interval = c(-0.1, 0), maximum = TRUE, alpha = alpha, eta = eta)$maximum)
}

#' Generate MxN matrix of SDE trajectories by monte carlo
#' Uses the Cox Ingersol Ross model:
#' dX[t] = (theta[1] - theta[2]*X[t])*dt + theta[3]sqrt(X[t])dW[t]
#' The CIR model is positive given the condition:
#' X[t] > 0 <=> 2*theta[1]*theta[2] >= theta[3]^2
#' @param param Parameters of a CIR model
#' @param M number of trajectories
#' @param N Number of data points per trajectory
#' @param X0 starting value of each trajectory
#' @param T.horizon Time length of each trajectory
#' @author Stav Sioutis
#' @rdname MC.CIR
#' @return MxN matrix
#' @export
MC.CIR = function(param, M, N, X0, T.horizon=N)
{
  # Matrix will be MxN
  cir = t(sde.sim(X0 = X0, theta=param, model="CIR", N=N, T=T.horizon, M=M))
  return(cir)
}

#' Fit Cox Ingersol Ross model:
#' dX[t] = (theta[1] - theta[2]*X[t])*dt + theta[3]sqrt(X[t])dW[t]
#' The CIR model is positive given the condition:
#' X[t] > 0 <=> 2*theta[1]*theta[2] >= theta[3]^2
#' @param series Numeric input series on which to fit the CIR model
#' @author Stav Sioutis
#' @rdname fit.CIR
#' @return parameters theta
#' @export
fit.CIR = function(series, verbose=FALSE)
{
  CIR.logistic = function(theta1, theta2, theta3)
  {
    if (!verbose)
      cat(theta1, "  ", theta2, "  ",theta3, "  \n")
    return(-sum(dsCIR(x=series, theta=c(theta1,theta2,theta3),log=TRUE)))
  }
  mu = mean(series)
  sd = sd(series)
  theta = mle(CIR.logistic, method='L-BFGS-B',
              start=list(theta1 = mu, theta2 = 0.1, theta3 = sd),
              lower=c(mu-sd, 0.001, 0.01),
              upper=c(mu+sd, 2, 1))
  return(coef(theta))
}

#' Generate MxN matrix of SDE trajectories by monte carlo
#' Uses the OU model:
#' dX[t] = (theta[1] - theta[2]*X[t])*dt + theta[3]dW[t]
#' @param param Parameters of the OU model
#' @param M number of trajectories
#' @param N Number of data points per trajectory
#' @param X0 starting value of each trajectory
#' @param T.horizon Time length of each trajectory
#' @author Stav Sioutis
#' @rdname MC.OU
#' @return MxN matrix
#' @export
MC.OU = function(param, M, N, X0, T.horizon=N)
{
  # Matrix will be MxN
  ou = t(sde.sim(X0 = X0, theta=param, model="OU", N=N, T=T.horizon, M=M))
  return(ou)
}

# Fit OU model
# dX[t] = (theta[1] - theta[2]*X[t])*dt + theta[3]dW[t]
#' @param series Numeric input series on which to fit the CIR model
#' @author Stav Sioutis
#' @rdname fit.CIR
#' @return parameters theta
#' @export
fit.OU = function(series)
{
  OU.logistic = function(theta1, theta2, theta3)
  {
    cat(theta1, "  ", theta2, "  ",theta3, "  \n")
    return(-sum(dsOU(x=series, theta=c(theta1,theta2,theta3),log=TRUE)))
  }
  mu = mean(series)
  sd = sd(series)
  theta = mle(OU.logistic, method='L-BFGS-B',
              start=list(theta1 = mu, theta2 = 0.1, theta3 = sd),
              lower=c(mu-sd, 0.001, 0.001),
              upper=c(mu+sd, 2, 1))
  return(coef(theta))
}

#' Generate MxN matrix of SDE trajectories by monte carlo
#' Uses the Double OU model:
#' S1[t] = lambda(mu - S1[t-1])dt + V1[t]*dW[t]
#' V1[t] = kappa(sigma - V1[t-1])dt + rho*dX[t]
#' @param param Parameters of the double OU model: c(lambda, kappa, rho)
#' @param mu Mean of the model
#' @param sigma Standard deviation of the model
#' @param M number of trajectories
#' @param N Number of data points per trajectory
#' @param X0 starting value of each trajectory
#' @param dt Time step
#' @author Stav Sioutis
#' @rdname MC.DoubleOU
#' @return MxN matrix
#' @export
MC.DoubleOU = function(param, mu, sigma, M=1000, N=500, X0, dt=1)
{
  lambda = param[1]
  kappa = param[2]
  rho = param[3]
  MC = matrix(0, M, N)
  for (i in 1:M)
  {
    # Simulate Double OU process
    dWt = sqrt(dt) * rnorm(N)
    dXt = sqrt(dt) * rnorm(N)
    S1 = rep(0, N)
    V1 = rep(0, N)
    S1[1] = mu
    V1[1] = rho
    for (t in 2:N)
    {
      V1[t] = V1[t-1] + kappa*(sigma-V1[t-1])*dt + rho*dXt[t]
      S1[t] = S1[t-1] + lambda*(mu-S1[t-1])*dt + V1[t]*dWt[t]
    }
    MC[i,] = S1
  }
  return(MC)
}

#' Fit Double OU process by minimizing histogram distance (chi-squared, weighted chi-squared, KL divergence, ...)
#' S1[t] = lambda(mu - S1[t-1])dt + V1[t]*dW[t]
#' V1[t] = kappa(sigma - V1[t-1])dt + rho*dX[t]
#' mu = mean(series) and sigma = sd(series) are fixed to improve the fitting procedure
#' @param series Numeric input series on which to fit the double OU model
#' @param M number of trajectories
#' @param N Number of data points per trajectory
#' @param mu mean mu of double OU model
#' @param sigma sd sigma of double OU model
#' @param lower lower bounds for each parameter
#' @param upper upper bounds for each parameter
#' @param chi.sq.w weighting parameter for weighted chi squared distance
#' @author Stav Sioutis
#' @rdname fit.DoubleOU
#' @return parameters c(lambda, kappa, rho)
#' @export
fit.DoubleOU = function(series, M=1000, N=300, 
                        mu = mean(series), sigma = sd(series),
                        lower=c(0.01, 0.01, 0.001),
                        upper=c(3, 1, 0.02),
                        chi.sq.w = 0.5)
{
  s.density = density(series)
  center = s.density$x[which.max(s.density$y)]
  
  obj = function(p)
  {
    ou = MC.DoubleOU(p, mu, sigma, tail(series,1), M, N)
    ou.density = density(ou)
    SS = 0
    cat(p[1], p[2], p[3], '\n')
    # Weighted chi sq (focuses on fitting tails)
    return(weighted.ChiSq(s.density, ou.density, chi.sq.w))
  }
  param = optim((upper-lower)/2,
                obj,
                method='L-BFGS-B',
                lower=lower,
                upper=upper)$par
  # S1[t] = lambda(mu - S1[t-1])dt + V1[t]*dW[t]
  # V1[t] = kappa(sigma - V1[t-1])dt + rho*dX[t]
  # param = c(lambda, kappa, rho)
  return(param)
}


#' Fit Double OU process by minimizing histogram distance (chi-squared, weighted chi-squared, KL divergence, ...)
#' S1[t] = lambda(mu - S1[t-1])dt + V1[t]*dW[t]
#' V1[t] = kappa(sigma - V1[t-1])dt + rho*dX[t]
#' mu = mean(series) and sigma = sd(series) are fixed to improve the fitting procedure
#' @param series Numeric input series on which to fit the double OU model
#' @param M number of trajectories
#' @param N Number of data points per trajectory
#' @param mu mean mu of double OU model
#' @param sigma sd sigma of double OU model
#' @param lower lower bounds for each parameter
#' @param upper upper bounds for each parameter
#' @param chi.sq.w weighting parameter for weighted chi squared distance
#' @author Stav Sioutis
#' @rdname fit.DoubleOU
#' @return parameters c(lambda, kappa, rho)
#' @export
fit.DoubleOU.par = function(series, M=1000, N=300, 
                        mu = mean(series), sigma = sd(series),
                        lower=c(0.01, 0.01, 0.001),
                        upper=c(3, 1, 0.02),
                        chi.sq.w = 0.5,
                        cores = min(4, detectCores()))
{
  obj = function(p)
  {
    ou = MC.DoubleOU(p, mu, sigma, tail(series,1), M, N)
    ou.density = density(ou)
    
    # Weighted chi sq (focuses on fitting tails)
    return(weighted.ChiSq(s.density, ou.density, chi.sq.w))
  }
  
  s.density = density(series)
  center = s.density$x[which.max(s.density$y)]
  
  cl = makeCluster(cores)
  registerDoParallel(cl)
  
  # split the first parameter space into separate regions for individual
  # cores to operate over
  bounds = seq(lower[1], upper[1], length.out = cores+1)
  
  # lower and upper bounds for each core
  low = matrix(0, nrow=cores, ncol=3)
  upp = matrix(0, nrow=cores, ncol=3)
  for (i in 1:cores)
  {
    low[i,] = c(bounds[i], lower[2], lower[3])
    upp[i,] = c(bounds[i+1], upper[2], upper[3])
  }
  
  param = foreach(i=1:cores, .combine = data.frame, 
                  .export=c('weighted.ChiSq', 'MC.DoubleOU')) %dopar% 
                  {
    j = optim((upp[i,]-low[i,])/2,
          obj,
          method='L-BFGS-B',
          lower=low[i,],
          upper=upp[i,],
          control = list('maxit' = floor(100/cores)))
    c(j$par, j$value)
                  }
  
  # S1[t] = lambda(mu - S1[t-1])dt + V1[t]*dW[t]
  # V1[t] = kappa(sigma - V1[t-1])dt + rho*dX[t]
  # param = c(lambda, kappa, rho)
  
  stopCluster(cl)
  return(param)
}


#' Calculate value of long position with txn.cost over a set of Monte Carlo trajectories
#' @param MC MxN input series of M Monte Carlo trajectories, each of length N
#' @param val Trade entry point
#' @param txn.cost Transaction cost
#' @param order.size Size of order to be placed. 
#' @author Stav Sioutis
#' @rdname long.pos
#' @return Expected return value of trade
long.pos = function(MC, val, txn.cost, order.size = 10)
{
  value = 0
  trades = 0
  for (i in 1:dim(MC)[1])
  {
    # First time the series crosses the long limit
    t1 = head(which(MC[i,] <= val), 1)
    t1 = ifelse(length(t1)==0, 0, t1)
    trades = trades + 1
    if (t1 == 0)
      next
    # Final value of the series
    final.val = tail(MC[i,],1)
    trade.val = final.val - MC[i,t1]
    value = value + trade.val * order.size - txn.cost
  }
  return(value/trades)
}

#' Calculate value of short position with txn.cost over a set of Monte Carlo trajectories
#' @param MC MxN input series of M Monte Carlo trajectories, each of length N
#' @param val Trade exit point
#' @param txn.cost Transaction cost
#' @param order.size Size of order to be placed. 
#' @author Stav Sioutis
#' @rdname short.pos
#' @return Expected return value of trade
short.pos = function(MC, val, txn.cost, order.size = 10)
{
  value = 0
  trades = 0
  for (i in 1:dim(MC)[1])
  {
    # First time the series crosses the long limit
    t1 = head(which(MC[i,] >= val), 1)
    t1 = ifelse(length(t1)==0, 0, t1)
    trades = trades + 1
    if (t1 == 0)
      next
    # Final value of the series
    final.val = tail(MC[i,],1)
    trade.val = MC[i,t1] - final.val
    value = value + trade.val * order.size - txn.cost
  }
  return(value/trades)
}

#' Find optimal entry/exit positions using matrix of MC trajectories
#' @param MC MxN input series of M Monte Carlo trajectories, each of length N
#' @param txn.cost Transaction cost
#' @param mu Upper bound of long positions and lower bound of short positions
#' @param sig Width of universe of long positions. Effectively the lower bound of long positions, and upper bound of short positions
#' @param seq.length Number of trade entry/exit points the MC algorithm with evaluate
#' @param order.size Size of order to be placed. 
#' @author Stav Sioutis
#' @rdname MC.eval
#' @return list of MC trajectories, along with optimal entry (long) and exit (short) positions
#' @export
MC.eval = function(MC, txn.cost=0.1, mu=mean(MC), sig=sd(MC)*3, seq.length=100, order.size=10)
{
  long = seq(mu, mu-sig, length.out = seq.length)
  best = 0
  best.idx = 0
  for (i in 1:length(long))
  {
    temp = long.pos(MC, long[i], txn.cost, order.size)
    if (temp > best)
    {
      best = temp
      best.idx = i
    }
  }
  long.opt = long[best.idx]
  
  short = seq(mu, mu+sig, length.out = seq.length)
  best = 0
  best.idx = 0
  for (i in 1:length(short))
  {
    temp = short.pos(MC, short[i], txn.cost, order.size)
    if (temp > best)
    {
      best = temp
      best.idx = i
    }
  }
  short.opt = short[best.idx]
  return(list("MC" = MC, "long" = long.opt, "short" = short.opt))
}

#' Plot estimated/emprical densities, cointegrating series, Monte Carlo simulations
#' and optimal entry/exit points
#' @param series cointegrated series vector
#' @param MC MxN input series of M Monte Carlo trajectories, each of length N
#' @param opt object returned by MC.eval containing entry and exit positions
#' @param num.trajectories number of MC trajectories to plot
#' @author Stav Sioutis
#' @rdname plot.MC
#' @export
plot.MC = function(series, MC, opt, num.trajectories = 50)
{
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1,2))
  par(mar=c(1.9,1.9,1,0.1))
  d.s = density(series)
  d.mc = density(MC)
  xlim=c(min(d.s$x, d.mc$x), max(d.s$x, d.mc$x))
  ylim=c(min(d.s$y, d.mc$y), max(d.s$y, d.mc$y))
  plot(d.s, main=paste('Distance:', weighted.ChiSq(d.s, d.mc)), xlim=xlim, ylim=ylim)
  lines(d.mc, col=2)
  legend('topleft', col=c(1,2), legend=c('Empirical', 'Fitted'), bty='n', lty=1)
  M = dim(MC)[1]
  N = dim(MC)[2]
  plot(series, type='l', main = paste('Cointegrating Series', q.stat(series)),
       xlim = c(0, length(series)+N),
       ylim = c(min(series, min(MC)), max(series, max(MC))))
  bleh = apply(MC[1:num.trajectories,], 1, FUN = lines, x=c(1:N)+length(series), col=rgb(0.5,0.5,0.5, 0.1))
  lines(t(MC)[,1], x=c(1:N)+length(series), col=2)
  abline(h = c(opt$short, opt$long), col=4, lty=3, lwd=2)
  legend('bottomright', legend=c('Max Expected Return'), col=c('blue'), lty=3, bty='n')
}









################## MONTE CARLO
# M = 1000
# S = coint
# vec = unlist(OUML(S, dt, N))
# mu = vec[1]
# sigma = vec[2]
# lambda = vec[3]
# 
# cat(sigma, lambda, mu)
# 
# MC = matrix(0, M, N)
# for (i in 1:M)
# {
#   # Simulate OU process
#   dWt = sqrt(dt) * rnorm(N)
#   S1 = rep(0, N)
#   S1[1] = tail(S,1)
# 
#   for (t in 2:N)
#     S1[t] = S1[t-1] + lambda*(mu-S1[t-1])*dt + sigma*dWt[t]
# 
#   MC[i,] = S1
#   lines(S1, x=(1:N)+length(S),col=rgb(0,0,1,0.05))
# }




