library(sn)

xi = 1
omega = 1
alpha = 1

nu = 1000
2*(m2+xi^2-2*xi*m1)/(m2+xi^2-2*xi*m1-omega^2)

delta = alpha/sqrt(1+alpha^2)
ifelse(nu>1e6,sqrt(pi/2),sqrt(pi/nu)*exp(lgamma(nu/2)-lgamma(nu/2-1/2)))*(m1-xi)/omega

mu = ifelse(nu<=1e6, delta*sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2)), delta*sqrt(2/pi) )
integrate( function(x){dst(x, xi, omega, alpha, nu)}, lower=-Inf, upper=Inf )
integrate( function(x){x*dst(x, xi, omega, alpha, nu)}, lower=-Inf, upper=Inf )

m1 = omega*mu+xi; m1
integrate( function(x){x^2*dst(x, xi, omega, alpha, nu)}, lower=-Inf, upper=Inf )

m2 = omega^2*ifelse(nu>1e6,1,nu/(nu-2))+2*xi*m1-xi^2; m2
integrate( function(x){x^3*dst(x, xi, omega, alpha, nu)}, lower=-Inf, upper=Inf )

m3 = omega^3*mu*(3-delta^2)*ifelse(nu>1e6,1,nu/(nu-3))+3*xi*m2-3*xi^2*m1+xi^3; m3
omega^2*(3-delta^2)*ifelse(nu>1e6,1,nu/(nu-3))*(m1-xi)+3*xi*m2-3*xi^2*m1+xi^3
omega^2*(m1-xi)*ifelse(nu>1e6,1,nu/(nu-3))*(3-(m1-xi)^2/omega^2*ifelse(nu>1e6,pi/2,pi/nu*exp(2*lgamma(nu/2)-2*lgamma(nu/2-1/2))))+3*xi*m2-3*xi^2*m1+xi^3
(m1-xi)*(3*omega^2*ifelse(nu>1e6,1,nu/(nu-3))-pi*ifelse(nu>1e6,1,nu/(nu-3))/2*(m1-xi)^2+3*pi/(4*(nu-3))*(m1-xi)^2)+3*xi*m2-3*xi^2*m1+xi^3
(m1-xi)*(3*omega^2*ifelse(nu>1e6,1,nu/(nu-3))-pi*ifelse(nu>1e6,1,nu/(nu-3))/2*(m1-xi)^2+3*pi/(4*(nu-3))*(m1-xi)^2)+3*xi*m2-3*xi^2*m1+xi^3
(m1-xi)*(m2+xi^2-2*xi*m1-omega^2)/(2*(m2-2*xi*m1+xi^2)-3*(m2+xi^2-2*xi*m1-omega^2))*(6*omega^2*(m2-2*xi*m1+xi^2)/(m2+xi^2-2*xi*m1-omega^2)-pi*((m2-2*xi*m1+xi^2)/(m2+xi^2-2*xi*m1-omega^2)-3/4)*(m1-xi)^2) +3*xi*m2-3*xi^2*m1+xi^3
(m1-xi)/(-m2+2*xi*m1-xi^2+3*omega^2)*(6*omega^2*(m2-2*xi*m1+xi^2)-pi*((m2-2*xi*m1+xi^2)-3/4*(m2+xi^2-2*xi*m1-omega^2))*(m1-xi)^2) +3*xi*m2-3*xi^2*m1+xi^3

integrate( function(x){x^4*dst(x, xi, omega, alpha, nu)}, lower=-Inf, upper=Inf )
m4 = 3*omega^4*ifelse(nu>1e6,1,nu^2/((nu-2)*(nu-4)))+4*xi*m3-6*xi^2*m2+4*xi^3*m1-xi^4; m4
3*omega^4*(2*m2+2*xi*m1+xi^2)/((2*m2+2*xi*m1+xi^2-2*(m2-omega^2))*(2*m2+2*xi*m1+xi^2-4*(m2-omega^2)))+4*xi*m3-6*xi^2*m2+4*xi^3*m1-xi^4
3*omega^4*(2*m2+2*xi*m1+xi^2)/((2*xi*m1+xi^2+2*omega^2)*(-2*m2+2*xi*m1+xi^2+4*omega^2))+4*xi*m3-6*xi^2*m2+4*xi^3*m1-xi^4
3*omega^2*(m2-2*xi*m1+xi^2)^2/(-m2+2*xi*m1-xi^2+2*omega^2)+4*xi*m3-6*xi^2*m2+4*xi^3*m1-xi^4


mom = function(m1, m2, m3, m4){
  omega2 = function(xi){
    C = (m4-4*xi*m3+6*xi^2*m2-4*xi^3*m1+xi^4)/(3*(m2-2*xi*m1+xi^2)^2)
    return( C*(-m2+2*xi*m1-xi^2)/(1-2*C) )
  }
  err_m3 = function(xi){
    C = m2 - 2*xi*m1 +xi^2
    m3Est = (m1-xi)/(-C+3*omega2(xi))*(6*omega2(xi)*C-pi*(C-3/4*(C-omega2(xi)))*(m1-xi)^2) +3*xi*m2-3*xi^2*m1+xi^3
    return( m3 - m3Est )
  }
  xiEst = uniroot( err_m3, interval=c(-10,10) )$root
  omegaEst = sqrt( omega2(xiEst) )
  nuEst = 2*(m2-2*xiEst*m1+xiEst^2) / (m2+xiEst^2-2*xiEst*m1-omegaEst^2)
  muEst = (m1-xiEst)/omegaEst
  deltaEst = muEst/ifelse(nuEst>1e6, sqrt(2/pi), sqrt(nu/pi)*exp(lgamma(nuEst/2-1/2)-lgamma(nuEst/2)) )
  alphaEst = deltaEst/sqrt(1-deltaEst^2)
  return( list(xi=xiEst, omega=omegaEst, alpha=alphaEst, nu=nuEst, mu=muEst, delta=deltaEst) )
}

xi = runif(1)
omega = rchisq(1,df=1)
alpha = runif(1)
nu = 10^runif(1,log(4)/log(10),3)
delta = alpha/sqrt(1+alpha^2)
mu = ifelse(nu<=1e6, delta*sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2)), delta*sqrt(2/pi) )
m1 = omega*mu+xi; m1
m2 = omega^2*ifelse(nu>1e6,1,nu/(nu-2))+2*xi*m1-xi^2; m2
m3 = omega^3*mu*(3-delta^2)*ifelse(nu>1e6,1,nu/(nu-3))+3*xi*m2-3*xi^2*m1+xi^3; m3
m4 = 3*omega^4*ifelse(nu>1e6,1,nu^2/((nu-2)*(nu-4)))+4*xi*m3-6*xi^2*m2+4*xi^3*m1-xi^4; m4
est = mom(m1, m2, m3, m4)
cat("xi   : ", round(xi,4), " ", round(est$xi,4), "\n", sep="")
cat("omega: ", round(omega,4), " ", round(est$omega,4), "\n", sep="")
cat("alpha: ", round(alpha,4), " ", round(est$alpha,4), "\n", sep="")
cat("nu   : ", round(nu,4), " ", round(est$nu,4), "\n", sep="")
