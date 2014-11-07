omega2 = function(xi){
  C = (m4-4*xi*m3+6*xi^2*m2-4*xi^3*m1+xi^4)/(3*(m2-2*xi*m1+xi^2)^2)
  return( C*(-m2+2*xi*m1-xi^2)/(1-2*C) )
}

err_m3 = function(xi){
  C = m2 - 2*xi*m1 +xi^2
  m3Est = (m1-xi)/(-C+3*omega2(xi))*(6*omega2(xi)*C-pi*(C-3/4*(C-omega2(xi)))*(m1-xi)^2) +3*xi*m2-3*xi^2*m1+xi^3
  return( m3 - m3Est )
}

#Function from package sn.
#I believe it takes four "moments" and returns the skew-t parameters.
#xbar = mean(d)
#s = sqrt(mean((d-xbar)^2))
#gamma1 <- mean((d-xbar)^3)/s^3
#gamma2 <- mean((d-xbar)^4)/s^4 - 3
#cp = c(xbar, s, gamma1, gamma2)
st.cp2dp <- function(cp, silent=FALSE, tol=1e-8, trace=FALSE) 
{
  fn0 <- function(log.nu, g1) st.gamma1(1, exp(log.nu)) - g1
  if(any(is.na(cp))) stop("NA's in argument 'cp'")
  p <- length(cp)-3
  x.names <- if(p>1) names(cp[2:p]) else NULL
  gamma1 <- cp[p+2]
  abs.g1 <- abs(gamma1)
  gamma2 <- cp[p+3]
  if(abs.g1 <=  0.5*(4-pi)*(2/(pi-2))^1.5)
    feasible <- (gamma2 > 2*(pi-3)*(2*abs.g1/(4-pi))^4/3)
  else {
    if(abs.g1 >= 4) feasible <- FALSE else {
      r0 <- uniroot(fn0, interval=c(log(4),1000), tol=tol, g1=abs.g1)
      nu0 <- exp(r0$root) 
      feasible <- (gamma2 >= st.gamma2(1,nu0))
      }
    }
  if(!feasible) {
    if(silent) return(NULL) else stop("CP outside feasible region")}
  delta <- 0.75*sign(gamma1)
  old <- c(delta,Inf)
  step <- Inf
  fn1 <- function(delta, g1, nu) st.gamma1(delta, nu) - g1
  fn2 <- function(log.nu, g2, delta) st.gamma2(delta, exp(log.nu)) - g2
  while(step > tol){
    fn21 <- fn2(log(4), gamma2, delta)
    fn22 <- fn2(log(100), gamma2, delta)
    if(any(is.na(c(fn21,fn22)))) stop("parameter inversion failed") # browser()
    if(fn21 * fn22 > 0)  return(rep(NA, p+3))
    r2 <- uniroot(fn2, interval=c(log(4),100), tol=tol, g2=gamma2, delta=delta)
    nu <- exp(r2$root)
    if(fn1(-1, gamma1, nu) * fn1(1, gamma1, nu)> 0) return(rep(NA, p+3))
    r1 <- uniroot(fn1, interval=c(-1,1), tol=tol, g1=gamma1, nu=nu)
    delta <- r1$root
    new <- c(delta, nu)
    step <- max(abs(old-new))
    if(trace) cat("delta, nu, log(step):", format(c(delta, nu, log(step))),"\n")
    old<- new
    }
  mu.z <- delta*b(nu)
  omega <- cp[p+1]/sqrt(nu/(nu-2) - mu.z^2)
  alpha <- delta/sqrt(1-delta^2)
  dp <- c(cp[1]-omega*mu.z, if(p>1) cp[2:p] else NULL, omega, alpha, nu)
  names(dp) <- param.names("DP", "ST", p, x.names=x.names)
  return(dp)
}

#My attempt at moments->parameters.  Doesn't seem to be stable.
mom = function(m1, m2, m3, m4){
  xiEst = uniroot( err_m3, interval=c(-10,10) )$root
  omegaEst = sqrt( omega2(xiEst) )
  #Force nu to be at least 5.  If nu=4, 4th moment DNE
  nuEst = min(5, 2*(m2-2*xiEst*m1+xiEst^2) / (m2+xiEst^2-2*xiEst*m1-omegaEst^2))
  muEst = (m1-xiEst)/omegaEst
  deltaEst = muEst/ifelse(nuEst>1e6, sqrt(2/pi), sqrt(nuEst/pi)*exp(lgamma(nuEst/2-1/2)-lgamma(nuEst/2)) )
  alphaEst = deltaEst/sqrt(1-deltaEst^2)
  return( list(xi=xiEst, omega=omegaEst, alpha=alphaEst, nu=nuEst, mu=muEst, delta=deltaEst) )
}

winsorizedEst = function(X, k=1.5){
  #k is traditionally sd's from the mean.  Convert to confidence interval assuming N(0,1).
  #Call it "confAlpha" since alpha is one of the params to estimate
  confAlpha = 1-pnorm(k)
  
  #Initial guess is N(median, mean absolute deviation):
  xi = median(X)
  omega = mad(X)
  alpha = 0
  nu = Inf
  
  converged = F
  while(!converged){
    lower = qst( confAlpha, xi, omega, alpha, nu )
    upper = qst( 1-confAlpha, xi, omega, alpha, nu )
    Y = X
    Y[Y>upper] = upper
    Y[Y<lower] = lower
    
    #Estimate moments
    ybar = mean(Y)
    res = Y-ybar
    s = sqrt(mean(res^2))
    gamma1 = mean(res^3)/s^3
    gamma2 = mean(res^4)/s^4-3
    fit = st.cp2dp( c(ybar, s, gamma1, gamma2), silent=TRUE, tol=1 )

    #With winsorized data, moments will be off.  Correct assuming the current dist.
    ybar = ybar*integrate( function(x){x*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val /
            integrate( function(x){min(max(x,lower),upper)*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val
    s = s*integrate( function(x){x^2*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val /
            integrate( function(x){min(max(x,lower),upper)^2*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val
    m3 = m3*integrate( function(x){x^3*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val /
            integrate( function(x){min(max(x,lower),upper)^3*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val
    m4 = m4*integrate( function(x){x^4*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val /
            integrate( function(x){min(max(x,lower),upper)^4*dst(x,xi,omega,alpha,nu)}, lower=-Inf, upper=Inf )$val
  
    est = mom(m1, m2, m3, m4)
    xi = est$xi
    omega = est$omega
    alpha = est$alpha
    nu = est$nu
    
    xGrid = seq(min(X), max(X), length.out=100)
    toPlot = dst(xGrid, xi, omega, alpha, nu)
    toPlot = data.frame( y=toPlot, x=xGrid )
    binWidth = diff(range(X))/20
    print( ggplot(data.frame(X)) + geom_bar(aes(x=X, weight=1/(binWidth*length(X))), binwidth=binWidth) +
      geom_line(data=toPlot, aes(x=x, y=y) ) )
    readline("Next")
  }
}

xiT = runif(1)
omegaT = rchisq(1,df=1)
alphaT = runif(1)
nuT = 10^runif(1,log(4)/log(10),3)
X = rst(1000, xi=xiT, omega=omegaT, alpha=alphaT, nu=nuT )
