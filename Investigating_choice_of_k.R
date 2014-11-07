#Currently, we're bounding the negative log-likelihood.  However, that's essentially a
#function of the density.  For the normal, the density goes down as sigma increases.  Thus,
#k will depend greatly on the true sigma.  Can we correct this?  Possibly choose a good value
#for k in the user inputs a hypothesized percent of outlier contamination? 

library(ggplot2)
library(plyr)
library(reshape)

#We need our robust estimator functions:
psi = function(x, k){
  if(x>k)
    return(exp(-x/k+1))
  return(1)
}

Psi = function(x, k){
#  if(k>=0)
#    ifelse(x>k, 2*k-k*exp(-x/k+1), x)
#  #If k<0, original function won't work.  Here's an alternative, but we might want
#  # something else...
#  else
#    ifelse(x>k,k*exp(x/k-1),x)
  ifelse(x<k, x, k+1-exp(-(x-k)) )
}

#Normal likelihood and robust likelihood functions
negLogLik = function(x, mu, sigma){
  1/2*log(2*pi) + log(sigma) + (x-mu)^2/(2*sigma^2)
}
negLogLikRob = function(x, mu, sigma, k){
  nLL = sapply(x, negLogLik, mu=mu, sigma=sigma)
  ifelse(nLL>20+k, k+1, Psi(nLL, k))
}
#Should be approximately invariant under scale
negLogLikRob2 = function(x, mu, sigma, k){
  nLL = 1/2*log(2*pi) + log(sigma) + Psi((x-mu)^2/(2*sigma^2), k=k)
  return(nLL)
}
#Attempt to reproduce the trimmed mean
negLogLikTrim = function(x, mu, sigma, alpha){
  z = (x-mu)/sigma
  k = -qnorm(alpha)
  nLL = 1/2*log(2*pi) + 1/2*log(sigma^2) + ifelse(z^2<k^2, z^2/2, k^2/2+sqrt(z^2/2-k^2/2))/2
  return(nLL)
}
#Attempt to reproduce the huber estimator
negLogLikWin = function(x, mu, sigma, k){
  z = (x-mu)/sigma
  nLL = 1/2*log(2*pi) + 1/2*log(sigma^2) + ifelse(z^2>k^2, k^2, z^2)/2 
  return(nLL)
}
gradNegLogLikWin = function(x, mu, sigma, k){
  z = (x-mu)/sigma
  dNLL = lapply(z, function(t){
    c(ifelse(t^2>k^2, 0, -2*t/sigma^2)/2 
    ,1/sigma + ifelse(t^2>k^2, 0, -2*t^2/sigma)/2)
  })
  dNLL = do.call("rbind", dNLL)
  dNLL = colSums(dNLL)
  return(dNLL)
}
x = rnorm(10)
mean(x); sd(x)
gradNegLogLikWin(x, mu=0, sigma=1, 10)
library(numDeriv)
grad(function(theta){sum(negLogLikWin(x, theta[1], theta[2], k=10))}, x=c(0,1) )

negLogLikRob(30, 0, 10, 3)/negLogLik(30, 0, 10)
negLogLikRob(3,  0,  1, 3)/negLogLik( 3, 0, 1)
negLogLikRob(.3, 0, .1, 3)/negLogLik(.3, 0, .1)
negLogLikRob2(30, 0, 10, 3)/negLogLik(30, 0, 10)
negLogLikRob2(3,  0,  1, 3)/negLogLik( 3, 0, 1)
negLogLikRob2(.3, 0, .1, 3)/negLogLik(.3, 0, .1)
#No winsorizing
negLogLikWin(30, 0, 10, 3)/negLogLik(30, 0, 10)
negLogLikWin(3,  0,  1, 3)/negLogLik( 3, 0, 1)
negLogLikWin(.3, 0, .1, 3)/negLogLik(.3, 0, .1)
#winsorizing
negLogLikWin(30, 0, 10, 2)/negLogLik(30, 0, 10)
negLogLikWin(3,  0,  1, 2)/negLogLik( 3, 0, 1)
negLogLikWin(.3, 0, .1, 2)/negLogLik(.3, 0, .1)

sigma = 10
xSeq = seq(-3,3,by=.01)*sigma
qplot( xSeq, dnorm(xSeq, sd=sigma) )

nSim = 100
kVals = c(seq(-5,10,.5),1000)
sigVals = c(.1,1,10,100,1000)

##########################################################################
# Initial simulation, goal is to show dependence of estimator on sigma
##########################################################################

estimates = matrix(ncol=4, nrow=0)
colnames(estimates)=c("sigma", "k", "mu_hat", "sigma_hat")
for(i in 1:nSim){
  x = rnorm(1000, mean=1)
  x[1:30] = rnorm(30, mean=6)
  
  for(sigma in sigVals){
    y = x*sigma
    temp = lapply( kVals, function(k){
      param = constrOptim(c(mean(y),sd(y)), f=function(theta){
        sum(negLogLikRob(y, theta[1], theta[2], k=k))
      }, ui=c(0,1), ci=0, grad=NULL )
      return( c(sigma=sigma, k=k,param$par) )
    } )
    temp = do.call("rbind", temp)
    estimates = rbind(estimates, temp)
  }
  cat("Simulation",i,"completed.\n")
}

estimates = data.frame(estimates)
ggplot( estimates, aes(x=k, y=mu_hat/sigma, color=as.factor(sigma), group=paste0(k,sigma)) ) +
  geom_boxplot() + labs(color="Sigma")

estimates$k = factor(estimates$k, levels=unique(c(estimates$k,"MSE")))
estimates$k[estimates$k==1000] = "MLE"
ggsave("Choice_of_k.png",
  ggplot( estimates, aes(x=k, y=mu_hat/sigma, color=as.factor(sigma), group=paste0(k,sigma)) ) +
    geom_boxplot() + labs(color="Sigma", y="mu_hat/mu")
)


##########################################################################
# Try fitting new negLogLikRob, one that penalizes only (x-mu)^2/sigma^2
##########################################################################

estimates2 = matrix(ncol=4, nrow=0)
colnames(estimates2)=c("sigma", "k", "mu_hat", "sigma_hat")
for(i in 1:nSim){
  x = rnorm(1000, mean=1)
  x[1:30] = rnorm(30, mean=6)
  
  for(sigma in sigVals){
    y = x*sigma
    temp = lapply( kVals, function(k){
      param = constrOptim(c(mean(y),sd(y)), f=function(theta){
        sum(negLogLikRob2(y, theta[1], theta[2], k=k))
      }, ui=c(0,1), ci=0, grad=NULL )
      return( c(sigma=sigma, k=k,param$par) )
    } )
    temp = do.call("rbind", temp)
    estimates2 = rbind(estimates2, temp)
  }
  cat("Simulation",i,"completed.\n")
}

estimates2 = data.frame(estimates2)
ggplot( estimates2, aes(x=k, y=mu_hat/sigma, color=as.factor(sigma), group=paste0(k,sigma)) ) +
  geom_boxplot() + labs(color="Sigma") +
  coord_cartesian(ylim=c(-.5,1.5))

estimates2$k = factor(estimates$k, levels=unique(c(estimates2$k,"MSE")))
estimates2$k[estimates$k==1000] = "MLE"
ggsave("Choice_of_k_partial_psi.png",
  ggplot( estimates2, aes(x=k, y=mu_hat/sigma, color=as.factor(sigma), group=paste0(k,sigma)) ) +
    geom_boxplot() + labs(color="Sigma", y="mu_hat/mu")
)

ggsave("Choice_of_k_partial_psi_zoomed.png",
  ggplot( estimates2, aes(x=k, y=mu_hat/sigma, color=as.factor(sigma), group=paste0(k,sigma)) ) +
    geom_boxplot() + labs(color="Sigma", y="mu_hat/mu") +
    coord_cartesian(ylim=c(0,2))
)

##########################################################################
# Compare my robust estimator with other common robust estimators
##########################################################################

d = data.frame( x = seq(-5,5,.1) )
k = 3
d$Squared = d$x^2/2
d$Absolute = abs(d$x)
d$Huber = ifelse(abs(d$x)<k, 1/2*d$x^2, k*(abs(d$x)-1/2*k) ) #http://en.wikipedia.org/wiki/Huber_loss
#Integral of Tukey biweight function
d$Tukey = ifelse(abs(d$x)<k, d$x^6/(6*k^4)-d$x^4/(2*k^2)+d$x^2/2, k^2/6)
d$New = Psi(d$x^2/2, k)
toPlot = melt(d, id.vars="x")
ggsave("~/Professional Files/Mines/Research/Robust Estimators/Comparison_of_M_estimators.png",
  ggplot( toPlot, aes(x=x, y=value, color=variable, group=variable)) + geom_line() +
    labs(x="x-hat{mu}", y="rho(x) (i.e. loss function)", color="Estimator", title="k=3")
  ,width=8, height=8)