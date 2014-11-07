library(sn)
library(ggplot2)

getPval = function(data, dp, probs=100:0/100){
  probs = sort(probs, decreasing=T)

  nu <- dp$nu
  Omega <- dp[[2]]
  Omega.bar <- cov2cor(Omega)
  alpha <- dp[[3]]
  alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
  omega <- sqrt(diag(Omega))
  l.nu <- (-1.3/nu - 4.93)
  h <- 100 * log(exp(((1.005 * alpha.star - 0.045) * l.nu - 
      1.5)/alpha.star) + 1)
  K <- h * (1.005 * alpha.star - 0.1) * (1 + nu)/(alpha.star * 
      nu)
  qF <- qf(probs, 2, nu)
  log.levels <- (lgamma(nu/2 + 1) - lgamma(nu/2) - log(pi * 
      nu) - 0.5 * log(1 - Omega.bar[1, 2]^2) - (nu/2 + 
      1) * log(2 * qF/nu + 1) + K - sum(log(omega)))

  dens = dmst(data, dp=dp)
  pval = rep(0, NROW(data))
  for(i in 1:length(probs)){
    filt = dens>exp(log.levels[i])
    pval[filt] = 1-probs[i]
  }
  return(pval)
}

out = rmst( 10000, dp=st2@dp )
pvalues = getPval( out, dp=st2@dp, probs=0:100/100 )
qplot( pvalues, binwidth=.01 )
sum( pvalues==1 )
