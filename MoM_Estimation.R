library(nleqslv)
library(BB)
library(ggplot2)

#Approximating gamma function ratio of nu:
nu = seq(4,1000,1)
qplot( nu, sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2)) ) +
  geom_line(aes(y=sqrt(nu/pi)*1/sqrt(nu/2-3/4)), color="red", linetype=2)
nu = seq(4,10,.01)
qplot( nu, sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2)) ) +
  geom_line(aes(y=sqrt(nu/pi)*1/sqrt(nu/2-3/4)), color="red", linetype=2)

stMoM = function(m1, m2, m3, m4){
  residual <- function(p) {
    xi = p[1]
    w = 100*exp(p[2])/(1+exp(p[2])) #Enforces omega to be positive
    nu = 1000*exp(p[3])/(1+exp(p[3])) #Enforces nu to be positive
    r = c(m2 - (w^2*nu/(nu - 2) + 2*xi*m1 - xi^2)
         ,m3 - ((m1-xi)*(3*w^2*nu-pi*(nu/2-3/4)*(m1-xi))/(nu-3) + 3*xi*m2 - 3*xi^2*m1 + xi^3)
         ,m4 - (w^4*3*nu^2/((nu - 2)*(nu - 4)) + 4*xi*m3 - 6*xi^2*m2 + 4*xi^3*m1 - xi^4))
    r
  }
  p0 <- c(0,1,10)
#  fit = nleqslv(p0, residual, control=list(maxit=1000))
  fit = BBsolve(p0, residual); fit$par

  xiEst = fit$par[1]
  omegaEst = exp(fit$par[2])
  nuEst = exp(fit$par[3])
  muEst = (m1-xiEst)/omegaEst
  deltaEst = muEst*(lgamma(nuEst/2)-lgamma(nuEst/2-1/2))*sqrt(pi/nuEst)
  alphaEst = deltaEst/sqrt(1-deltaEst^2)
  return( c(xi=xiEst, omega=omegaEst, alpha=alphaEst, nu=nuEst, muEst, deltaEst) )
}

stMoM2 = function(m1, m2, m3, m4){
  xi = 0
  omega = 1
  alpha = 0
  nu = 10
  delta = alpha/sqrt(1+alpha^2)
  mu = delta*(lgamma(nu/2-1/2)-lgamma(nu/2))*sqrt(nu/pi)

  xi = m1 - omega*mu
  omega = sqrt( max(0.001,(nu-2)/nu*( m2-2*xi*m1+xi^2 )) )
  resid = function(delta){(m3-3*xi*m2+3*xi^2*m1-xi^3)/omega^3- delta*(lgamma(nu/2-1/2)-lgamma(nu/2))*sqrt(nu/pi)*(3 - delta^2)*nu/(nu - 3)}
  alpha = nleqslv(alpha, resid)$x
  resid = function(nu){m4 - 3*omega^4*nu^2/((nu - 2)*(nu - 4)) + 4*xi*m3 - 6*xi^2*m2 + 4*xi^3*m1 - xi^4}
  nu = max(4,nleqslv(nu, resid)$x)
  delta = alpha/sqrt(1+alpha^2)
  mu = delta*(lgamma(nu/2-1/2)-lgamma(nu/2))*sqrt(nu/pi)
}

n = 10000
xiSam = runif(n,-10,10)
omegaSam = 10^runif(n,-2,2)
alphaSam = runif(n,-10,10)
nuSam = 10^runif(n,log(4)/log(10),4)
delta = alphaSam/sqrt(1+alphaSam^2)
mu = delta*(lgamma(nuSam/2-1/2)-lgamma(nuSam/2))*sqrt(nuSam/pi)

m1 = omegaSam*mu + xiSam
m2 = omegaSam^2*nuSam/(nuSam - 2) + 2*xiSam*m1 - xiSam^2
m3 = omegaSam^3*mu*(3 - delta^2)*nuSam/(nuSam - 3) + 3*xiSam*m2 - 3*xiSam^2*m1 + xiSam^3
m4 = 3*omegaSam^4*nuSam^2/((nuSam - 2)*(nuSam - 4)) + 4*xiSam*m3 - 6*xiSam^2*m2 + 4*xiSam^3*m1 - xiSam^4

stMoM(m1, m2, m3, m4)

library(nnet)
y=cbind(xiSam, omegaSam, alphaSam, nuSam)
fit = nnet( x=cbind(m1,m2,m3,m4), y, size=20, linout=T, maxit=1000)
fit = nnet( x=cbind(m1,m2,m3,m4), y=xiSam, size=20, linout=T, maxit=100000 )
summary( fit$fitted.values - xiSam )



library(sn)
xi = 3
omega = 1
alpha = -1
nu = Inf
d = rst(1000, xi, omega, alpha, nu)
mean(d); sd(d)*length(d)/(length(d)-1)
stMoM(d)

stMoM = function(data){
  #Initial guesses: N(0,1)
  xiFit = 0
  omegaFit = 1
  alphaFit = 0
  nuFit = Inf

  converged = FALSE
  while(!converged){
    delta = alphaFit/sqrt(alphaFit^2+1)
    mu = delta*ifelse(nuFit>1000, sqrt(2/pi), sqrt(nuFit/pi)*exp(lgamma(nuFit/2-1/2)-lgamma(nuFit/2)) )

    #Fit via method of moments.  First estimate xi
    EY = mean(data)
    xiFit = EY-omegaFit*mu
    
    #Compute centered second moment and estimate omega
    EY2 = mean((data-xiFit)^2)
    omegaFit = EY2*ifelse(nuFit>1000,1,(nuFit-2)/nuFit)
  
    #Estimate alpha and nu using 3rd and 4th central, standardized moments
    EY3 = mean( ((data-xiFit)/omegaFit)^3 )
    f_of_nu = ifelse(nuFit>1000, sqrt(2/pi), nuFit/(nuFit-3)*sqrt(nuFit/pi)*gamma(nuFit/2-1/2)/gamma(nuFit/2) )
    #k = ( -EY3/(2*f_of_nu) + 1/2*sqrt(EY3^2/f_of_nu^2-4))^(1/3)
    deltaFit = nleqslv(0, function(x){x*(3-x^2)-EY3/f_of_nu})$x
    alphaFit = deltaFit/sqrt(1-deltaFit^2)
    EY4 = mean( ((data-xiFit)/omegaFit)^4 )
    nuFit = ifelse(EY4<3, Inf, (3*EY4+sqrt(EY4*(EY4+24)))/(EY4-3) ) #Negative solution is always bounded above by 2, but nu>4 so that doesn't work
    cat( round( c(xiFit, omegaFit, alphaFit, nuFit), 4 ), "\n" )
    Sys.sleep(0.2)
  }
  
  return(c(xiFit, omegaFit, alphaFit, nuFit))
}