#' Quantify stationarity
#' 
#' Quantify stationarity of a numeric input series
#' @param x input series
#' @param k1 First window width for calculating rolling mean/standard deviations
#' @param k2 Second window width for calculating rolling mean/standard deviations
#' @author Stav Sioutis
#' @rdname q.stat
#' @export
q.stat = function(x, k1=floor(0.5*length(x)), k2=floor(0.25*length(x)))
{
  mu = mean(x)
  sig = sd(x)
  # mean zero variance 1
  z = (x - mu)/sig
  # sum squared differences of moving average and moving sd
  sma = SMA(z, n = k1) + SMA(z, n=k2)
  ssd = SMA(z^2, n = k1) + SMA(z^2, n=k2) - 2
  stat = sum(na.omit(sma^2 + ssd^2))
  return(stat/sqrt(length(x)))
}



#' Find optimal entry/exit positions using matrix of MC trajectories
#' @param density1 First density object returned by \code{density}
#' @param density2 Second density object returned by \code{density}
#' @param weight Weight applied to the tails (higher weight encourages tail similarity)
#' @author Stav Sioutis
#' @rdname weighted.Chisq
#' @return distance between two densities
#' @export
weighted.ChiSq = function(density1, density2, weight = 0.5)
{
  SS = 0
  chisq = 0
  n = length(density1$y)
  center = density1$x[which.max(density1$y)]
  for (i in 1:n)
  {
    chisq = ((density1$y[i] - density2$y[i])^2)/((density1$y[i] + density2$y[i])^2)
    SS = SS + chisq * (weight*(density1$x[i]-center)^2+1)
  }
  return(SS/(2*n))
}



#' Load ETF from IBrokers
#' Load ETF along with holdings, construct and return NAV along with weights vector and matrix of holdings
#' @param etf ETF of interest
#' @param holdings csv file containing holdings and weights
#' @param duration length of data to retrieve
#' @param barSize frequency of data
#' @param wd directory containing csv file with holdings data
#' @author Stav Sioutis
#' @rdname load.ETF
#' @return ETF data 
#' @export
load.ETF = function(etf, holdings, duration = '1 Y', barSize = '1 day',
                    wd = 'C:/Users/Stavros/Google Drive/Trading/R/cointegration/')
{
  etf.holdings = read.csv(paste(wd, holdings, sep = ''), as.is=T, header = T, skip=1)
  etf.sym = sub('.', " ", toupper(etf.holdings[,1]), fixed = T)
  etf.wts = as.numeric(sub('%', '', etf.holdings[,3]))/100
  
  # Load every stock in the sector on which the ETF is defined
  etf.h = xts()
  etf.vol = xts()
  tws = twsConnect()
  count = 0
  for (i in etf.sym)
  {
    count = count+1
    cat(count,'/',length(etf.sym), ' ')
    data = reqHistoricalData(tws,
                             twsEquity(i, 'SMART', 'ISLAND'),
                             duration = duration, barSize=barSize)
    etf.h = merge(etf.h, data[,4])
  }
  
  data = reqHistoricalData(tws,
                          twsEquity(etf, 'SMART', 'ISLAND'),
                          duration = duration, barSize=barSize)
  ETF = data[,4]
  etf.vol = data[,5]
  twsDisconnect(tws)
  etf.h = na.omit(etf.h)
  pc = list("h" = etf.h, "vol" = etf.vol, "wts" = etf.wts, "etf" = tail(na.omit(ETF), dim(etf.h)[1]))
  return(pc)
}

#' Load Equities from IBrokers
#' @param sym Equities of interest
#' @param duration length of data to retrieve
#' @param barSize frequency of data
#' @author Stav Sioutis
#' @rdname load.Equity
#' @return Equity data
#' @export
load.Equity = function(sym, duration, bar)
{
  data = xts()
  tws = twsConnect()
  for (i in sym)
  {
    data = merge(data, reqHistoricalData(tws,
                                         twsEquity(i, 'SMART', 'ISLAND'),
                                         duration = duration, barSize = bar)[,4])
  }
  twsDisconnect(tws)
  data = na.omit(data)
  return(data)
}