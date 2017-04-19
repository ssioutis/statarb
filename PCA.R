library(IBrokers)
library(quantmod)


# source('util.R')

#' Principle components analysis of sector ETF
#' @param etf ETF
#' @param pca.trunc proportion of smallest principle components to drop
#' @param offset.idx remove first n observations
#' @author Stav Sioutis
#' @rdname pc.analysis
#' @export
pc.analysis = function(etf, pca.trunc = 0.1, offset.idx = 0)
{
  # Std dev grid lines (entry/exit points)
  sd.vec = c(-1.2, -0.5, 0.5, 1.2)
  
  len = dim(etf$etf)[1]
  
  etf$wts = tail(etf$wts, len - offset.idx)
  etf$h   = tail(etf$h, len - offset.idx)
  etf$etf = tail(etf$etf, len - offset.idx)
  
  
  cat(length(etf$wts), 'assets\n')
  # PCA
  pc = princomp(etf$h)
  # Plot eigenvalues
  screeplot(pc)
  # Loading vectors
  pc1 = pc$loadings[,1]
  pc2 = pc$loadings[,2]
  pc3 = pc$loadings[,3]
  # Remove insignificant contributions
  # pc1 = ifelse(abs(pc1) < pca.trunc, 0, pc1)
  # pc2 = ifelse(abs(pc2) < pca.trunc, 0, pc2)
  # pc3 = ifelse(abs(pc3) < pca.trunc, 0, pc3)
  # Each 'zero' means we are not trading a particular asset, more zeroes = less txn costs
  cat('PC1 zeros', length(which(pc1==0))/length(pc1), '\n')
  # Principal components
  pc1.vec = etf$h %*% pc1
  pc2.vec = etf$h %*% pc2
  pc3.vec = etf$h %*% pc3
  # Regress principal components on ETF
  pc1.m = lm(etf$etf ~ pc1.vec)
  pc2.m = lm(etf$etf ~ pc2.vec)
  pc3.m = lm(etf$etf ~ pc3.vec)
  # Scale PCs accordingly
  pc1.scale = coef(pc1.m)[1] + coef(pc1.m)[2]*pc1.vec
  pc2.scale = coef(pc2.m)[1] + coef(pc2.m)[2]*pc2.vec
  pc3.scale = coef(pc3.m)[1] + coef(pc3.m)[2]*pc3.vec
  
  # Factor model
  # factor = pc1.vec + pc2.vec + pc3.vec
  # factor.m = lm(etf$etf ~ factor)
  # factor.scale = coef(factor.m)[1] + coef(factor.m)[2]*factor
  
  ##################################### plot 2
  par(mfrow=c(2,1), mar=c(0,4,0,0))
  #Plot ETF with top 3 principal components
  plot(as.numeric(etf$etf), type='l', lwd=2, col=1, ylab='ETF/PCs', xaxt='n')
  lines(pc1.scale, col=2, lwd=2)
  lines(pc2.scale, col=3, lwd=2)
  lines(pc3.scale, col=4, lwd=2)
  legend('bottomleft', col=c(2,3,4,1), legend=c('pc1', 'pc2', 'pc3', 'ETF'), bty='n', lty=1, lwd=2)
  #Plot ETF/pc1 spread
  spread = as.numeric(etf$etf-pc1.scale)
  par(mar=c(2,4,0,0))
  plot(spread, type='l', lwd=2, ylab = 'ETF/PC spread')
  abline(h = mean(spread), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread) + mean(spread), col='grey', lty=2)
  
}

#' NAV analysis of sector ETF
#' @param etf ETF
#' @param force.drop drops the assets with the smallest weights first before computing q.stat
#' @author Stav Sioutis
#' @rdname NAV.analysis
#' @return Top 4 stationary spreads
#' @export
NAV.analysis = function(etf, force.drop = floor(0.5 * length(etf$wts)))
{
  # Std dev grid lines (entry/exit points)
  sd.vec = c(-1.2, -0.5, 0.5, 1.2)
  
  # Need to find truncation parameter which:
  # a) maximizes stationarity of the cointegrating series
  # b) maximizes expected return of an optimal OU strategy
  
  
  # Find NAV truncation parameter which maximizes stationarity of the spread
  n = length(etf$wts)
  # each row is a copy of etf$wts
  # first column is for q.stats
  # best.wts = data.frame(0, t(matrix(etf$wts, nrow=n, ncol=n-3)))
  best.wts = data.frame(0, t(matrix(etf$wts, nrow=n, ncol=n-force.drop-3)))
  m = dim(best.wts)[2]
  for (i in 1:(n-force.drop-3))
  {
    # zero out the next smallest weight
    best.wts[i,(m-i-force.drop):m] = 0
    # try new set of weights
    curr.wts = as.numeric(best.wts[i, 2:m])
    # recreate spread
    nav.t = etf$h %*% curr.wts
    # linear model
    nav.m.t = lm(as.numeric(etf$etf) ~ nav.t)
    # scaled NAV to fit over ETF
    nav.scale.t = coef(nav.m.t)[1] + coef(nav.m.t)[2] * nav.t
    # spread between scaled NAV and ETF
    spread.nav.t = as.numeric(etf$etf-nav.scale.t)
    q = q.stat(spread.nav.t)
    best.wts[i,1] = q
  }
  best.wts = best.wts[order(best.wts[,1]),]
  wts.t = as.numeric(best.wts[1,2:m])
  print(wts.t)
  
  # wts.t = ifelse(abs(etf$wts < nav.trunc), 0, etf$wts)
  
  cat('NAV zeros', length(which(wts.t==0))/length(wts.t), '\n')
  
  ##################################### plot 1
  # Plot ETF, NAV and NAV-truncated
  par(mfrow=c(3,1), mar=c(0,4,0,0))
  plot(as.numeric(etf$etf), type='l', lwd=2, ylab='ETF/NAV', xaxt='n')
  nav = etf$h %*% etf$wts
  
  nav.t = etf$h %*% wts.t
  nav.m = lm(as.numeric(etf$etf) ~ nav)
  nav.m.t = lm(as.numeric(etf$etf) ~ nav.t)
  nav.scale = coef(nav.m)[1] + coef(nav.m)[2] * nav
  nav.scale.t = coef(nav.m.t)[1] + coef(nav.m.t)[2] * nav.t
  lines(nav.scale, col=2, lwd=2)
  lines(nav.scale.t, col=3, lwd=2)
  legend('bottomleft', col=c(1,2,3), lwd=2, bty='n', legend=c('ETF', 'NAV', 'NAV-truncated'))
  # Plot ETF/NAV spread
  spread.nav = as.numeric(etf$etf-nav.scale)
  plot(spread.nav, type='l', lwd=2, ylab = 'ETF/NAV spread', xaxt='n')
  abline(h = mean(spread.nav), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread.nav) + mean(spread.nav), col='grey', lty=2)
  # Plot ETF volume in log scale
  par(mar=c(2,4,0,0))
  barplot(etf$vol, col=c(2,3), ylab = 'ETF Volume')
  
  ##################################### plot 2
  # Plot top 4 stationary truncated series 
  par(mfrow=c(4,1), mar=c(0,4,0,0))
  nav.t1 = etf$h %*% as.numeric(best.wts[1,2:m])
  nav.t2 = etf$h %*% as.numeric(best.wts[2,2:m])
  nav.t3 = etf$h %*% as.numeric(best.wts[3,2:m])
  nav.t4 = etf$h %*% as.numeric(best.wts[4,2:m])
  nav.m.t1 = lm(as.numeric(etf$etf) ~ nav.t1)
  nav.m.t2 = lm(as.numeric(etf$etf) ~ nav.t2)
  nav.m.t3 = lm(as.numeric(etf$etf) ~ nav.t3)
  nav.m.t4 = lm(as.numeric(etf$etf) ~ nav.t4)
  # Works better when we fit with intercept but exclude intercept in the
  # reconstruction of the scaled series
  nav.scale.t1 = coef(nav.m.t1)[2] * nav.t1
  nav.scale.t2 = coef(nav.m.t2)[2] * nav.t2
  nav.scale.t3 = coef(nav.m.t3)[2] * nav.t3
  nav.scale.t4 = coef(nav.m.t4)[2] * nav.t4
  spread.nav.t1 = as.numeric(etf$etf-nav.scale.t1)
  spread.nav.t2 = as.numeric(etf$etf-nav.scale.t2)
  spread.nav.t3 = as.numeric(etf$etf-nav.scale.t3)
  spread.nav.t4 = as.numeric(etf$etf-nav.scale.t4)
  
  plot(spread.nav.t1, type='l', lwd=2, ylab = 'ETF/NAV-trunc spread', xaxt='n')
  legend('topleft', legend=paste('q.stat =', q.stat(spread.nav.t1)))
  abline(h = mean(spread.nav.t1), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread.nav.t1) + mean(spread.nav.t1), col='grey', lty=2)
  
  plot(spread.nav.t2, type='l', lwd=2, ylab = 'ETF/NAV-trunc spread', xaxt='n')
  legend('topleft', legend=paste('q.stat =', q.stat(spread.nav.t2)))
  abline(h = mean(spread.nav.t2), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread.nav.t2) + mean(spread.nav.t2), col='grey', lty=2)
  
  plot(spread.nav.t3, type='l', lwd=2, ylab = 'ETF/NAV-trunc spread', xaxt='n')
  legend('topleft', legend=paste('q.stat =', q.stat(spread.nav.t3)))
  abline(h = mean(spread.nav.t3), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread.nav.t3) + mean(spread.nav.t3), col='grey', lty=2)
  
  par(mar=c(2,4,0,0))
  plot(spread.nav.t4, type='l', lwd=2, ylab = 'ETF/NAV-trunc spread')
  legend('topleft', legend=paste('q.stat =', q.stat(spread.nav.t4)))
  abline(h = mean(spread.nav.t4), col='red', lty=2, lwd=2)
  abline(h = sd.vec * sd(spread.nav.t4) + mean(spread.nav.t4), col='grey', lty=2)
  
  
  #reset graphics
  par(mfrow=c(1,1))
  
  s1 = list('spread' = spread.nav.t1, 'NAV' = nav.m.t1)
  s2 = list('spread' = spread.nav.t2, 'NAV' = nav.m.t2)
  s3 = list('spread' = spread.nav.t3, 'NAV' = nav.m.t3)
  s4 = list('spread' = spread.nav.t4, 'NAV' = nav.m.t4)
  
  return(list('s1'=s1,'s2'=s2,'s3'=s3,'s4'=s4))
}




