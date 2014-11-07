#setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
source("Code/sn-funct.R")
source("Code/modified_mst.pdev.grad.R")

#Original mst.pdev.grad is not additive for D matrix.  New implementation fixes this.
mst.pdev.grad = mst.pdev.grad.new

#Function implementing the robustness adjustment to the density.  If x is some supplied
#   negative log-likelihood value, then psi(x) returns the adjusted negative log-likelihood.
#   Essentially, if the negative log-likelihood value is too large, it is bounded via this
#   function.
#x: negative log-likelihood
#k: parameter controlling the amount of robustification.  Negative log-likelihoods larger
#   than k are reduced, and the adjusted value is always less than 2*k.
psi = function(x, k){
  if(x>k)
    return(exp(-x/k+1))
  return(1)
}

#Function implementing the robustness adjustment to the derivative of the density.  If x is
#   some supplied negative log-likelihood value, then Psi(x) returns the derivative of the
#   adjusted negative log-likelihood.  Essentially, if the negative log-likelihood value is
#   too large, it is bounded via psi().  This function corrects the derivative.
#x: negative log-likelihood
#k: parameter controlling the amount of robustification.  Negative log-likelihoods larger
#   than k are reduced, and the adjusted value is always less than 2*k.
Psi = function(x, k){
  if(x>k)
    return(2*k-k*exp(-x/k+1))
  return(x)
}

#Let's try a simple case first: N(mu, sigma), both parameters unknown
n.dev = function(dp, y){
  mu = dp[1]
  sigma = dp[2]
  return( sum( (y-mu)^2/(2*sigma^2) ) + length(y)*log(sigma) )
}

#Deviance is "additive", i.e. we can add deviances together for total deviance
n.dev.robust = function(dp, y, k=2){
  nonRobust = lapply(y, n.dev, dp=dp)
  nonRobust = do.call("c", nonRobust)
  robust = ifelse(nonRobust>k, sapply(nonRobust, Psi, k=k), nonRobust)
  return( sum(robust) )
}

n.dev.gh = function(dp, y){
  mu = dp[1]
  sigma = dp[2]
  return( c( -sum(y-mu)/(sigma^2), -sum((y-mu)^2/sigma^3) + length(y)/sigma ) )
}

#Gradient is "additive", i.e. we can compute it for each obs and then add all estimates together:
n.dev.gh.robust = function(dp, y, k=2){
  nonRobustGH = lapply(y, n.dev.gh, dp=dp)
  nonRobustGH = do.call("rbind", nonRobustGH)
  nonRobust = lapply(y, n.dev, dp=dp)
  nonRobust = do.call("c", nonRobust)
  #If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't change for that case.
  #If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
  robustGH = nonRobustGH
  if(any(nonRobust>k))
    robustGH[nonRobust>k,] = robustGH[nonRobust>k,]*sapply(nonRobust[nonRobust>k], psi, k=k)
  return( colSums(robustGH) )
}

#Computes the robust-ified, penalized deviance for the univariate skew-t.
#dp: "Density parameters"?  Vector of xi, omega, alpha, nu.
#x: matrix of the independent variables for the fit.  Typically just a matrix of ones.
#y: vector of dependent variables
#k: Parameter controlling the robustness of the fit.  The largest possible value for the
#   negative log-likelihood is 2*k, and the negative log-likelihood is adjusted down 
#   whenever it is larger than k.
st.pdev.robust = function(dp, x, y, k=2, ...){
  nonRobust = lapply(y, st.pdev, dp=dp, x=matrix(1,nrow=1))
  nonRobust = do.call("c", nonRobust)
  robust = ifelse(nonRobust>k, sapply(nonRobust, Psi, k=k), nonRobust)
  return( sum(robust) )
}

#Computes the robust-ified, penalized deviance for the multivariate skew-t.
#param: Optimization parameters derived from xi, Omega, alpha, nu.  For conversion, see
#   the functions optpar2dplist and dplist2optpar.
#x: matrix of the independent variables for the fit.  Typically just a matrix of ones.
#y: matrix of dependent variables
#k: Parameter controlling the robustness of the fit.  The largest possible value for the
#   negative log-likelihood is 2*k, and the negative log-likelihood is adjusted down 
#   whenever it is larger than k.
mst.pdev.robust = function(param, x, y, k=2, ...){
  nonRobust = sapply(1:NROW(y), function(i){
    mst.pdev( param, x=matrix(1), y=y[i,,drop=F], w=1 )
  })
  robust = ifelse(nonRobust>k, sapply(nonRobust, Psi, k=k), nonRobust)
  return( sum(robust) )
}

#Computes the gradient of the robust-ified, penalized deviance for the univariate skew-t.
#dp: "Density parameters"?  Vector of xi, omega, alpha, nu.
#x: matrix of the independent variables for the fit.  Typically just a matrix of ones.
#y: vector of dependent variables
#k: Parameter controlling the robustness of the fit.  The largest possible value for the
#   negative log-likelihood is 2*k, and the negative log-likelihood is adjusted down 
#   whenever it is larger than k.
st.pdev.gh.robust = function(dp, x, y, k=2){
  nonRobustGH = lapply(y, st.pdev.gh, dp=dp, x=matrix(1,nrow=1))
  nonRobustGH = do.call("rbind", nonRobustGH)
  nonRobust = lapply(y, st.pdev, dp=dp, x=matrix(1,nrow=1))
  nonRobust = do.call("c", nonRobust)
  #If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't change for that case.
  #If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
  robustGH = nonRobustGH
  if(any(nonRobust>k))
    robustGH[nonRobust>k,] = robustGH[nonRobust>k,]*sapply(nonRobust[nonRobust>k], psi, k=k)
  #Some NA's can occur from st.pdev.gh.  Set them to 0, since I have no better idea as to why they occur (presumably really small density values)
  robustGH[is.na(robustGH)] = 0
  return( colSums(robustGH) )
}

#Computes the gradient of the robust-ified, penalized deviance for the multivariate skew-t (with respect
#   to the optimization parameters).
#param: Optimization parameters derived from xi, Omega, alpha, nu.  For conversion, see
#   the functions optpar2dplist and dplist2optpar.
#x: matrix of the independent variables for the fit.  Typically just a matrix of ones.
#y: matrix of dependent variables
#k: Parameter controlling the robustness of the fit.  The largest possible value for the
#   negative log-likelihood is 2*k, and the negative log-likelihood is adjusted down 
#   whenever it is larger than k.
mst.pdev.grad.robust = function(param, x, y, k=2, ...){
  nonRobustGH = t( sapply(1:NROW(y), function(i){
    mst.pdev.grad(param=param, x=x[i,,drop=F], y=y[i,,drop=F], w=1)
    }) )
  nonRobust = sapply(1:NROW(y), function(i){
    mst.pdev(param=param, x=x[i,,drop=F], y=y[i,,drop=F], w=1)
    })
  #If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't change for that case.
  #If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
  robustGH = nonRobustGH
  if(any(nonRobust>k))
    robustGH[nonRobust>k,] = robustGH[nonRobust>k,]*sapply(nonRobust[nonRobust>k], psi, k=k)
  #Some NA's can occur from st.pdev.gh.  Set them to 0, since I have no better idea as to why they occur (presumably really small density values)
  robustGH[is.na(robustGH)] = 0
  return( colSums(robustGH) )  
}

#n: Sample size to simulate
#outPct: Percent of outliers.  Actual number is round(n*outPct)
#outSigma: Outliers are created by adding on an error of N(0,outSigma^2)
#k: A numeric vector of constants to use in "capping" the likelihood function.
#  Each constant is tried in turn, and results are returned for all values.
#xi, omega, alpha, nu: The skew-t parameters to simulate.
#pressure: Use Denver station data at this press. level, and simulate the skew-t using Ying's code
#type: The type of skewness to use.  Must be "MVN", "obs", or "EX"
#restrict: Should outliers be generated in the tail of the skew-t only?  Only applies if
#  pressure and type are not null.
runSim = function(n=100, outPct=0, outSigma=0, k=6:10, fast.k=c(.99,.98,.95,.9)
    ,xi0 = runif(1, min=-10, max=10)
    ,Omega0 = diag(x=10^runif(1,-1,1), nrow=1)
    ,alpha0 = runif(1,-4,4)
    ,nu0 = 10^runif(1,log(4)/log(10),4)
    ,pressure = NULL, type = NULL
    ,restrict=FALSE){ 

  d = length(xi0)
  p = 1 #p = ncol(x)

  #Input parameter checks:  
  if(!is.numeric(n))
    stop("n must be numeric!")
  if(n<10)
    stop("n must be >= 10 for even a chance at reasonable results!")
  if(outPct>1 | outPct<0 )
    stop("outPct must be between 0 and 1!")
  if(outSigma<0)
    stop("outSigma must be non-negative!")
  if(any(k<1))
    stop("Values of k less than 1 are likely to give very poor results.")
  if(any(dim(Omega0)!=p))
    stop("Omega0 must have dimension pxp, where p=length(xi0)!")
  if(length(alpha0)!=p)
    stop("alpha0 must have the same length as xi0!")
  if(length(nu0)!=1)
    stop("nu0 must have length 1!")
  if(!is.null(pressure))
    if(!pressure %in% c(700,500,400,300,250,200,100,70))
      stop("Invalid value for pressure!")
  if(!is.null(type))
    if(!type %in% c("obs", "MVN", "EX") )
      stop("Invalid value for type!")
  if(is.null(pressure) & !is.null(type))
    stop("pressure and type must both be specified or both be NULL!")
  if(!is.null(pressure) & is.null(type))
    stop("pressure and type must both be specified or both be NULL!")
  if(restrict & is.null(pressure))
    stop("restrict=TRUE is only implemented for non-NULL pressure and type arguments!")
  
  if(is.null(pressure)){
    #Generate random data:
    x = matrix(rep(1,n))
    p = length(xi0)
    if(p==1){ #Univariate data
      y = as.vector(rst(n, xi0, Omega0, alpha0, nu0))    
      #Contaminate with outliers?
      outRows = c()
      if(outPct>0){
        outRows = sample(n, size=round(outPct*n), replace=F)
        y[outRows] = ( y[outRows]
          + matrix(sample(c(-1,1),size=length(outRows),replace=T),ncol=1)
          * matrix(rnorm(length(outRows), mean=20, sd=outSigma), ncol=1) )
      }
    } else { #Multivariate data
      y = rmst(n, xi=xi0, Omega=Omega0, alpha=alpha0, nu=nu0)
      y = matrix( as.numeric(y), ncol=2 )
      #Contaminate with outliers?
      outRows = c()
      if(outPct>0){
        outRows = sample(n, size=round(outPct*n), replace=F)
        y[outRows,] = ( y[outRows,]
          + matrix(sample(c(-1,1),size=2*length(outRows),replace=T),ncol=2)
          * matrix(rnorm(2*length(outRows), mean=20, sd=outSigma), ncol=2) )
      }
    }
  } else { #pressure is specified, generate data using Ying's code.
    denver = denverParams(type=type)
    r = (1:16)[grepl(pressure,names(denver$xi))]
    params = marginal( denver$xi, denver$omega, denver$alpha, denver$nu, r )
    #Put parameters in workspace so we can save to simulation results
    xi0 = params$xi; Omega0=params$omega; alpha0=params$alpha; nu0=params$nu
    #Note: the k below is not the same k as in the function arguments!  That k is a tuning
    #parameter for the robust fitting algorithm.
    if(!restrict)
      data = m2Angle(n, p=1, xi=xi0, omega=Omega0, alpha=alpha0, nu=nu0, p.out=outPct, k=runif(n,10,12) )
    else
      data = m2AngleRestrict(n, p=1, xi=xi0, omega=Omega0, alpha=alpha0, nu=nu0, p.out=outPct, k=runif(n,10,12) )

    #Construct bivariate skew-t data for fitting
    x = matrix(1, nrow=n)
    y = cbind(data$u, data$v)
    outRows = (1:NROW(x))[data$outlier]
  }
  
  d = NCOL(y)
  p = ncol(x)
  
  rnames = NULL
  
  start = Sys.time()
  nlmEst = robustST(y, robust=F, method="constrOptim", w=rep(1,nrow(x)), k=k[1])
  if(d>1) #Multivariate fit returns a list, need to concatenate to a vector
    nlmEst = c(nlmEst$beta, as.numeric(nlmEst$Omega), nlmEst$alpha, nlmEst$nu, nlmEst$convergence )
  #Bind on the computation time
  nlmEst = c(nlmEst, as.numeric(difftime( Sys.time(), start, units="mins")) )
  nlmEst = c(nlmEst, NA)
  est = nlmEst
  rnames = c(rnames, "nlmEst")
  
  start = Sys.time()
  conEst = robustST(y, robust=F, method="constrOptim", w=rep(1,nrow(x)), k=k[1])
  if(d>1) #Multivariate fit returns a list, need to concatenate to a vector
    conEst = c(conEst$beta, as.numeric(conEst$Omega), conEst$alpha, conEst$nu, conEst$convergence )
  #Bind on the computation time
  conEst = c(conEst, as.numeric(difftime( Sys.time(), start, units="mins")) )
  conEst = c(conEst,NA)
  est = rbind(est, conEst)
  rnames = c(rnames, "conEst")

  for(kVal in k){
    start = Sys.time()
    nlmRobEst = robustST(y, robust=T, method="nlminb", w=rep(1,nrow(x)), k=kVal)
    if(d>1) #Multivariate fit returns a list, need to concatenate to a vector
      nlmRobEst = c(nlmRobEst$beta, as.numeric(nlmRobEst$Omega), nlmRobEst$alpha, nlmRobEst$nu, nlmRobEst$convergence )
    #Bind on the computation time
    nlmRobEst = c(nlmRobEst, as.numeric(difftime( Sys.time(), start, units="mins")) )
    nlmRobEst = c(nlmRobEst, kVal)
    est = rbind(est, nlmRobEst)
    rnames = c(rnames, "nlmRobEst")
    
    start = Sys.time()
    conRobEst = robustST(y, robust=T, method="constrOptim", w=rep(1,nrow(x)), k=kVal)
    if(d>1) #Multivariate fit returns a list, need to concatenate to a vector
      conRobEst = c(conRobEst$beta, as.numeric(conRobEst$Omega), conRobEst$alpha, conRobEst$nu, conRobEst$convergence )
    #Bind on the computation time
    conRobEst = c(conRobEst, as.numeric(difftime( Sys.time(), start, units="mins")) )
    conRobEst = c(conRobEst, kVal)
    est = rbind(est, conRobEst)
    rnames = c(rnames, "conRobEst")
  }
  
  if(NCOL(y) %in% 1:2){
    if(NCOL(y)==1)
      optFunc = fast.TLE.ST
    else
      optFunc = fast.TLE.BST
    
    for(kVal in fast.k){
      start = Sys.time()
      TLEest = try( optFunc(y, k=kVal) )
      if(d>1 & !is(TLEest,"try-error")) #Multivariate fit returns a list, need to concatenate to a vector
        TLEest = c(TLEest$xi, as.numeric(TLEest$omega), TLEest$alpha, TLEest$nu, NA )
      if(d==1 & !is(TLEest,"try-error")) #Multivariate fit returns a list, need to concatenate to a vector
        TLEest = TLEest$MLE
      if(is(TLEest,"try-error"))
        TLEest = rep(NA, ifelse(d==1,4,9)) #univariate has 4 parameters, bivariate has 9
      #Bind on the computation time
      TLEest = c(TLEest, as.numeric(difftime( Sys.time(), start, units="mins")) )
      TLEest = c(TLEest, kVal)
      est = rbind(est, TLEest)
      rnames = c(rnames, "TLEest")
    }
  }
      
  #Update column names for xi, omega, alpha, nu.  Loops are used since p is a variable
  cnames =c()
  for(i in 1:d)
    cnames = c(cnames, paste0("xi_",i))
  for(i in 1:d)
    for(j in 1:d)
      cnames = c(cnames, paste0("omega_",i,j))
  for(i in 1:d)
    cnames = c(cnames, paste0("alpha_",i))
  cnames =c(cnames, "nu", "convergence", "compTime", "k")
  colnames(est) = cnames
  
  est = data.frame(est)
  est$modelName = rnames
  est$n = n
  est$outCnt = length(outRows)
  #Add in columns for true xi, omega, alpha, nu.  Loops are used since p is a variable
  for(i in 1:d)
    est = cbind(est, xi0[i])
  for(i in 1:d)
    for(j in 1:d)
      if(d==1){
        est = cbind(est, Omega0)
      } else {
        est = cbind(est, Omega0[i,j])
      }
  for(i in 1:d)
    est = cbind(est, alpha0[i])
  est = cbind(est, nu0)
  
  #Set the column names for the true parameters to the column names for the estimated ones.
  colnames(est)[(d^2+2*d+1+6):(2*d^2+4*d+2+6-1)] = paste0(colnames(est)[1:(d^2+2*d+1)],"_act")

  return(est)
}

#Fits a robust version of the multivariate skew-t, done by bounding the negative log-likelihood
#  for each observation.
#y: vector or matrix of observations to fit the skew-t to.
#x: matrix of ones, or matrix of independent variables for skew-t regression (currently untested)
#robust: Should the robust estimator be used?
#method: constrOptim uses a constrained algorithm, forcing nu and omega>0.  However, the
#  implementation for multivariate skew-t fitting enforces this by default, so nlminb and
#  constrOptim should be very similar for multivariate data.  For univariate, constrOptim
#  is recommended.  For multivariate, constrOptim is also recommended as it appears to be
#  faster.
#w: vector of case weights, defaults to 1.
#k: parameter controlling the "robustness" of the fit.  The maximum value for the negative
#  log-likelihood for any observation is 2*k.  Thus, as k->Inf the estimator approaches the
#  MLE.  k values around 8 or 10 seem to perform well.
#start: starting values for the optimization
robustST = function(y, x=matrix(1,nrow=NROW(y)), robust=T, method=c("nlminb", "constrOptim"), w=rep(1,nrow(x)), k=10
            ,start=NULL){
  library(sn)
  
  #Data quality checks
  if(any(is.na(y))){
    if(is.null(dim(y)))
      filt = !is.na(y)
    else
      filt = !apply(y, 1, function(x){any(is.na(x))})
    
    x = x[filt,]
    w = w[filt]
    if(is.null(dim(y)))
      y = y[filt]
    else
      y = y[filt,]
  }
  if(!is(x,"matrix"))
    stop("x must be a matrix!")
  if(!is.matrix(y) & !is.numeric(y))
    stop("y must be a matrix or numeric vector!")
  if(!is(robust,"logical"))
    stop("robust must be a logical!")
  if(nrow(x)!=NROW(y))
    stop("x and y must have the same number of observations!")
  if(!is.numeric(w))
    stop("w must be numeric!")
  if(length(w)!=nrow(x))
    stop("w must have the same length as ncol(x)!")
  if(length(method)>1)
    method = method[1]
  if(!method %in% c("nlminb", "constrOptim"))
    stop("method must be one of nlminb or constrOptim!")
  
  n = nrow(x)
  p = ncol(x)
  d = NCOL(y)
  nw = sum(w)
  
  #Treat univariate separately
  if(d==1){
    if(is.null(start)){
      #Determine starting estimate (via logic from st.mple function in sn)
      ls <- lm.wfit(x, y, w)
      res <- ls$residuals
      s <- sqrt(sum(w * res^2)/nw)
      gamma1 <- sum(w * res^3)/(nw * s^3)
      gamma2 <- sum(res^4)/(nw * s^4) - 3
      cp <- c(ls$coef, s, gamma1, gamma2)
      dp <- st.cp2dp(cp, silent = TRUE)
      if (is.null(dp)) 
          dp <- rep(NA, length(cp))
      if (any(is.na(dp))) 
          dp <- c(cp[1:(p + 1)], 0, 10)
      names(dp) = c("xi", "omega", "alpha", "nu")
    } else {
      dp = start
    }
    
    if(!robust & method=="nlminb"){
      nlmEst = try(nlminb( start=dp
              ,function(dp){st.pdev(dp, x, y, w=w)}
              ,gradient=function(dp){st.pdev.gh(dp, x, y)}))
      if(is(nlmEst, "try-error")){
        fit = rep(NA, length(dp)+1)
      } else {
        fit = c( nlmEst$par, convergence=nlmEst$convergence )
      }
    }
    
    if(!robust & method=="constrOptim"){
      conEst = try(constrOptim(theta=dp
              ,f=function(dp){st.pdev(dp, x, y)}
              ,grad=function(dp){st.pdev.gh(dp, x, y)}
              #Constraints: force omega>0 and nu>0
              ,ui=matrix(c(0,0,1,0,0,0,0,1), nrow=2)
              ,ci=rep(0,2)))
      if(is(conEst, "try-error")){
        fit = rep(NA, length(dp)+1)
      } else {
        fit = c( conEst$par, convergence=conEst$convergence )
      }
    }
  
    if(robust & method=="nlminb"){
      nlmRobEst = try(nlminb( start=dp
              ,function(dp){st.pdev.robust(dp, x, y, k=k)}
              ,gradient=function(dp){st.pdev.gh.robust(dp, x, y, k=k)}))
      if(is(nlmRobEst, "try-error")){
        fit = rep(NA, length(dp)+1)
      } else {
        fit = c( nlmRobEst$par, convergence=nlmRobEst$convergence )
      }
    }
    
    if(robust & method=="constrOptim"){
      conRobEst = try(constrOptim(theta=dp
              ,f=function(dp){st.pdev.robust(dp, x, y, k=k)}
              ,grad=function(dp){st.pdev.gh.robust(dp, x, y, k=k)}
              ,ui=matrix(c(0,0,1,0,0,0,0,1), nrow=2)
              ,ci=rep(0,2)))
      if(is(conRobEst, "try-error")){
        fit = rep(NA, length(dp)+1)
      } else {
        fit = c( conRobEst$par, convergence=conRobEst$convergence )
      }
    }
    
    return(fit)
    
  } else { #Now multivariate case
    if(is.null(start)){
      #Determine starting estimate (via logic from mst.mple function in sn)
      ls <- lm.wfit(x, y, w, singular.ok=FALSE)
      beta <- coef(ls)
  		Omega <-  var(resid(ls))
  		omega <- sqrt(diag(Omega))
  		alpha <- rep(0, d)
      nu <- 8
      param <- dplist2optpar(list(beta=beta, Omega=Omega, alpha=alpha))
      param <- c(param, log(nu))
    } else {
      param = start
    }
  
    if(!robust & method=="nlminb"){
      nlmEst = try(nlminb( start=param
              ,function(param){mst.pdev(param, x, y)}
              ,gradient=function(dp){mst.pdev.grad(param, x, y, w=rep(1,NROW(x)))}))
      if(is(nlmEst, "try-error")){
        fit = rep(NA, length(param)+1)
      } else {
        fit = c(nlmEst$par, nlmEst$convergence)
      }
    }
    
    if(!robust & method=="constrOptim"){
      conEst = try(constrOptim(theta=param
            ,f=function(param){mst.pdev(param, x, y, w=w)}
            ,grad=function(param){mst.pdev.grad(param, x, y, w=w)}
            #No need for constraints as optpar2dplist ensures nu>0 and Omega is pos. def.
            #So, set u_i to all 0's, and force this to always be greater than -1
            # (which it always will).
            ,ui=matrix(0,ncol=length(param))
            ,ci=-1))
      if(is(conEst, "try-error")){
        fit = rep(NA, length(param)+1)
      } else {
        fit = c(conEst$par, conEst$convergence)
      }
    }
  
    if(robust & method=="nlminb"){
      nlmRobEst = try(nlminb( start=param
              ,function(param){mst.pdev.robust(param, x, y, k=k)}
              ,gradient=function(param){mst.pdev.grad.robust(param, x, y, k=k)}))
      if(is(nlmRobEst, "try-error")){
        fit = rep(NA, length(param)+1)
      } else {
        fit = c( nlmRobEst$par, nlmRobEst$convergence )
      }
    }
    
    if(robust & method=="constrOptim"){
      conRobEst = try(constrOptim(theta=param
              ,f=function(param){mst.pdev.robust(param, x, y, k=k)}
              ,grad=function(param){mst.pdev.grad.robust(param, x, y, k=k)}
              #No need for constraints as optpar2dplist ensures nu>0 and Omega is pos. def.
              #So, set u_i to all 0's, and force this to always be greater than -1
              # (which it always will).
              ,ui=matrix(0,ncol=length(param))
              ,ci=-1))
      if(is(conRobEst, "try-error")){
        fit = rep(NA, length(param)+1)
      } else {
        fit = c(conRobEst$par, conRobEst$convergence)
      }
    }
    
    optpar = optpar2dplist(fit[-length(fit)], p=p, d=d)
    return( list(beta=optpar$beta, Omega=optpar$Omega, alpha=optpar$alpha, nu=optpar$nu, convergence=fit[length(fit)] ) )
  }
}

plot_results = function(results, prefix, xi=0, omega=1, alpha=0, nu=10000){
  toPlot = results
  toPlot$n = paste0("n=", toPlot$n)
  #Ensure ordering is appropriate
  toPlot$n = factor(toPlot$n, levels=paste0("n=", sort(unique(results$n))) )
  toPlot$outPct = paste0("Outliers: ",toPlot$outPct,"%")
  #Ensure ordering is appropriate
  toPlot$outPct = factor(toPlot$outPct, levels=paste0("Outliers: ", sort(unique(results$outPct)), "%") )
  
  myTheme = theme( axis.text.x=element_text(size=14)
                  ,axis.title.x=element_text(size=14)
                  ,axis.text.y=element_text(size=14)
                  ,axis.title.y=element_text(size=14)
                  ,strip.text.x=element_text(size=14)
                  ,strip.text.y=element_text(size=14))
  
  p = ggplot(toPlot, aes(y=xiEst, x=estimator) ) + geom_boxplot() +
      facet_grid( n ~ outPct ) + geom_hline(yintercept=xi, color="red", linetype=2) +
      myTheme
  ggsave( paste0(prefix, "_xi_est.png"), p, width=8, height=10 )
  
  p = ggplot(toPlot, aes(y=omegaEst, x=estimator) ) + geom_boxplot() +
      facet_grid( n ~ outPct ) + geom_hline(yintercept=omega, color="red", linetype=2) +
      myTheme
  ggsave( paste0(prefix, "_omega_est.png"), p, width=8, height=10 )

  p = ggplot(toPlot, aes(y=alphaEst, x=estimator) ) + geom_boxplot() +
      facet_grid( n ~ outPct ) + geom_hline(yintercept=alpha, color="red", linetype=2) +
      myTheme
  ggsave( paste0(prefix, "_alpha_est.png"), p, width=8, height=10 )

  p = ggplot(toPlot, aes(y=nuEst, x=estimator) ) + geom_boxplot() +
      facet_grid( n ~ outPct ) + geom_hline(yintercept=nu, color="red", linetype=2) +
      scale_y_log10(breaks=10^(1:10)) + myTheme
  ggsave( paste0(prefix, "_nu_est.png"), p, width=8, height=10 )
}

#Convert multivariate parameters into bivariate using formulas from wind_radiosonde_QC_01.pdf (from Mandy)
#alphaAdj: If TRUE, computes the correct alpha for the marginal (per Azzalini et.al's
# "Distributions generated by perturbations of symmetry", page 18 and Capitano et.al's
# "Graphical models for skew-normal variates", page 15).  If false, just subsets alpha.
marginal = function( xi, omega, alpha, nu, r=1, alphaAdj=F ){
  #Check data types
  if(!is.numeric(xi))
    stop("xi must be numeric!")
  if(!is(omega,"matrix"))
    stop("omega must be a matrix!")
  if(!is.numeric(alpha))
    stop("alpha must be numeric!")
  if(!is.numeric(nu))
    stop("nu must be numeric!")

  #Check data dimensions
  p = length(xi)
  s = 1:p
  s = s[-r]
  if(any(dim(omega)!=p))
    stop("omega must have dimensions pxp, where p=length(xi)!")
  if(length(alpha)!=p)
    stop("alpha must have the same length as xi")
  if(length(nu)!=1)
    stop("nu must have length 1!")
  if(any(!r %in% 1:p))
    stop("All of r's values must be in 1:p!")
  
  xi.r = xi[r]
  omega.r = omega[r,r]
  if(alphaAdj){
    w = diag(sqrt(diag(omega)))
    omega_b = diag(1/diag(w)) %*% omega %*% diag(1/diag(w)) #Note: diag(1/diag(w)) is w^{-1}
    w.r = diag(sqrt(diag(omega.r)))
    B = diag(1/diag(w))%*%omega[,r]
    alpha.r = w.r %*% solve(omega.r) %*% t(B) %*% alpha /
      as.numeric( sqrt(1+t(alpha)%*%(omega_b - B %*% solve(omega.r) %*% t(B)) %*% alpha) )
  } else {
    alpha.r = alpha[r]
  }
  nu.r = nu
  return( list(xi=xi.r, omega=omega.r, alpha=alpha.r, nu=nu.r) )
}

#This function approximates the integral of the abs. value of the difference in 2d skew-t densities.
diffEst = function(param, param2, xGridPts=30, yGridPts=xGridPts){
  grid = defineGrid(param, xGridPts, yGridPts)
  int = integrate2d(grid[[1]], grid[[2]], function(x){
    abs(dmst(x, xi=param[[1]], Omega=param[[2]], alpha=param[[3]], nu=param[[4]]) -
      dmst(x, xi=param2[[1]], Omega=param2[[2]], alpha=param2[[3]], nu=param2[[4]]))
  })
  return(int)
}

#This function chooses a grid for evaluating the 2d integral of the difference in densities.
defineGrid = function(param, xGridPts=30, yGridPts=xGridPts){
  xi = param[[1]]
  omega = param[[2]]
  alpha = param[[3]]
  nu = param[[4]]

  #Naive approach: assume bivariate skew-t is best represented by two marginal skew-t's.
  #xGrid = qst( 1:(xGridPts-1)/xGridPts, xi[1], omega[1,1], alpha[1], nu )
  #yGrid = qst( 1:(yGridPts-1)/yGridPts, xi[2], omega[2,2], alpha[2], nu )
  #return(list(xGrid, yGrid))
  
  #Naive approach 2: Step by omega and go xGridPts/2*omega or yGridPts/2*omega in +/- directions
  xGrid = param[[1]][1] + (-xGridPts/2):(xGridPts/2)*param[[2]][1,1]/2
  yGrid = param[[1]][2] + (-yGridPts/2):(yGridPts/2)*param[[2]][2,2]/2
  return(list(xGrid, yGrid))
}

#Integrates func() over x and y space.  The integral is estimated using a 2-d analogue of the
# trapesoidal rule.  Each "trapesoid" has a rectangular base and 4 (possibly different)
# heights of the corners.  We compute the volume as A(z1+z2+z3+z4)/4, where A is the base area
# and zi the ith height.  Note that trapesoids are only defined within the grid, thus all
# volumne outside of the grid is ignored.
#xGrid: the x values to approximate the integral over.
#yGrid: the y values to approximate the integral over.
#func: a function whose first argument is a 2d vector, and the function evaluates to a scalar.
#...: additional arguments to pass to func.
integrate2d = function( xGrid, yGrid, func, ...){
  library(reshape)
  
  grid = merge( xGrid, yGrid)
  colnames(grid) = c("x", "y")

  vals = apply( grid, 1, func, ... )
  data = melt( cbind(grid, vals), id.vars=1:3 )
  data = cast( data, x ~ y, value="vals")
  #Remove xGrid from data (just want data values in the matrix)
  rownames(data) = data[,1]
  data[,1] = NULL
  
  #Average columns of data
  data2 = matrix(0, nrow=nrow(data), ncol=nrow(data)-1)
  for( i in 1:(ncol(data)-1) ){
    data2[,i] = (data[,i] + data[,i+1])/2
  }
  
  data = data2
  data2 = matrix(0, nrow(data)-1, ncol=ncol(data) )
  for( i in 1:(nrow(data)-1) ){
    data2[i,] = (data[i,] + data[i+1,])/2
  }
  
  widths = diff( xGrid )
  heights = diff( yGrid )
  areas = t(t(widths))%*%t(heights)
  volumes = areas*data2
  return(sum(volumes))
}

#From ying_model3.R (Denver station parameters)
denverParams = function(type=c("MVN", "obs", "EX")){
  if(length(type)>1)
    type==type[1]
  if(!type %in% c("MVN", "obs", "EX"))
    stop("type must be one of 'MVN', 'obs', or 'EX'!")
  
  xi=c(3.5732607, 10.3322993, 13.0043508, 16.5135646, 19.0336571, 20.8985299, 11.4608094, 5.2606013, -1.0291708, -2.2748794, -1.7431219, -1.2281764, -0.9888074, -0.4555666, -0.7184110, -0.9835646)
  names(xi)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  
  omega=rbind( c(22.3972699, 13.984818,  16.039642,  17.684241,  18.63175,  18.346261, 10.4795136,  8.2677082, -3.9072388, -10.004767, -10.995055, -13.833277, -13.634031, -10.2153270, -2.815235, -0.5339653), c(13.9848178, 53.727990,  59.020103,  66.056794,  66.75493,  59.537442, 23.6810410, 14.8809412,  5.5295733,   3.477148,   6.406637,   8.474633,   9.738715,  11.3172025,  4.378073,  3.3395517), c(16.0396421, 59.020103,  90.528398, 102.709659, 103.20949,  88.314971, 30.7581432, 17.6257511,  9.3983732,  11.193745,  17.330145,  24.295520,  26.207537,  25.1415773,  8.529959,  6.3410444), c(17.6842407, 66.056794, 102.709659, 147.899721, 148.89594, 124.512607, 38.0446733, 19.9519976, 12.4206720,  18.832824,  27.023667,  38.891596,  42.174570,  38.6743609, 12.975820,  8.9548713), c(18.6317507, 66.754935, 103.209486, 148.895942, 171.79385, 145.080702, 43.3933979, 22.0252564, 12.3625644,  19.464066,  27.533566,  40.567982,  44.611103,  40.2599930, 13.941600,  9.4114303), c(18.3462610, 59.537442,  88.314971, 124.512607, 145.08070, 149.877677, 46.6103745, 24.3935571, 10.9785169,  14.172450,  21.046724,  30.593732,  35.229002,  33.2482152, 10.321343,  7.7821568), c(10.4795136, 23.681041,  30.758143,  38.044673,  43.39340,  46.610374, 51.6696619, 33.6064021,  0.1524421,  -4.419451,  -4.078705,  -3.998690,  -3.327270,  -0.0476168, -1.178374,  1.0679616), c( 8.2677082, 14.880941,  17.625751,  19.951998,  22.02526,  24.393557, 33.6064021, 40.3398453, -2.7411404,  -7.184680,  -7.805223,  -9.335120,  -9.578368,  -5.4294290, -2.678479, -0.6314295), c(-3.9072388,  5.529573,   9.398373,  12.420672,  12.36256,  10.978517,  0.1524421, -2.7411404, 25.9163894,  27.043738,  30.793393,  37.083542,  38.886431,  35.2738880, 13.268837,  6.7957864),c(-10.0047669,  3.477148,  11.193745,  18.832824,  19.46407,  14.172450, -4.4194513, -7.1846796, 27.0437381,  78.522421,  90.324931, 108.070579, 112.710362, 100.0875035, 36.757158, 17.5988472),c(-10.9950547,  6.406637,  17.330145,  27.023667,  27.53357,  21.046724, -4.0787046, -7.8052225, 30.7933930,  90.324931, 126.802112, 151.111722, 155.933320, 135.8397389, 47.631396, 22.8484928),c(-13.8332773,  8.474633,  24.295520,  38.891596,  40.56798,  30.593732, -3.9986897, -9.3351204, 37.0835417, 108.070579, 151.111722, 209.976779, 215.636087, 185.3566163, 62.552556, 29.9704567),c(-13.6340314,  9.738715,  26.207537,  42.174570,  44.61110,  35.229002, -3.3272702, -9.5783683, 38.8864311, 112.710362, 155.933320, 215.636087, 242.728801, 208.9461148, 70.680936, 33.6640122),c(-10.2153270, 11.317202,  25.141577,  38.674361,  40.25999,  33.248215, -0.0476168, -5.4294290, 35.2738880, 100.087504, 135.839739, 185.356616, 208.946115, 208.5533958, 72.232734, 34.7282020), c(-2.8152352,  4.378073,   8.529959,  12.975820,  13.94160,  10.321343, -1.1783744, -2.6784786, 13.2688371,  36.757158,  47.631396,  62.552556,  70.680936,  72.2327344, 46.452689, 21.1557030), c(-0.5339653,  3.339552,   6.341044,   8.954871,   9.41143,   7.782157,  1.0679616, -0.6314295,  6.7957864,  17.598847,  22.848493,  29.970457,  33.664012,  34.7282020, 21.155703, 18.8168900))
  
  colnames(omega)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  rownames(omega)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  
  #Constructing an omega that is perfectly symmetric.  omega above fails the symmetry test because of round-off error.
  UT=upper.tri(omega)
  omega2=matrix(0,nrow=nrow(omega),ncol=ncol(omega))
  omega2[UT]=omega[UT]
  omega2=omega2+t(omega2)
  diag(omega2)=diag(omega)
  eigen(omega2)$values
  colnames(omega2) = colnames(omega)
  rownames(omega2) = rownames(omega)
  
  #Values based on those observed at Denver Station
  alpha.obs=c(2.16,  1.44,  1.35,  0.95,  0.60,  0.31,  2.62,  3.02, -0.03, -0.81, -0.82, -0.89, -0.82, -0.75, -0.31,  0.39 )
  #Rewrite so we have a 1,-1 for the 300 pressure levels
  alpha.obs=c(2.16,  1.44,  1.35,  1,  0.60,  0.31,  2.62,  3.02, -0.03, -0.81, -0.82, -1, -0.82, -0.75, -0.31,  0.39 )
  #alpha.obs=c(0.09, -0.07,  0.13, -0.25,  0.18, -0.08, -0.04,  0.00, -0.09,  0.11,  0.02, -0.08, -0.03,  0.07,  0.07, -0.10)
  df.obs=10
  names(alpha.obs) = names(xi)

  #Values for a MVN distribution
  alpha.MVN=rep(0,16)
  df.MVN=Inf
  names(alpha.MVN) = names(xi)
  
  #Values for distributions that are more skewed than what was observed at the Denver Station
  alpha.EX=c(3:10,3:10)
  #Rewrite so we have a 6,0 for the 300 pressure levels
  alpha.EX[4]=6; alpha.EX[12]=0
  df.EX=5
  names(alpha.EX) = names(xi)
  
  if(type=="MVN")
    return(list( xi=xi, omega=omega2, alpha=alpha.MVN, nu=df.MVN) )
  if(type=="obs")
    return(list( xi=xi, omega=omega2, alpha=alpha.obs, nu=df.obs) )
  if(type=="EX")
    return(list( xi=xi, omega=omega2, alpha=alpha.EX, nu=df.EX) )
}








#############################################################################################
#
#
# Data simulation functions from Ying
#
#
#############################################################################################

#model 1 without outliers
m1=function(n,p,xi,omega,alpha,nu){
  #p by n data matrix
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  list(u=u,v=v)
}

#outlying at all levels
m2=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  C=rbinom(2*n,1,p.out)
  
  s=2*rbinom(2*n,1,0.5)-1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
m2Adj=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,2*n) #JB
  C[1:(p.out*n)] = 1 #JB
  C = sample(C) #JB
  
  s=2*rbinom(2*n,1,0.5)-1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  k = sqrt(rep( diag(omega), each=n))*k #JB
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away.  This function differs from m2Adj in that random
# errors (outliers) are generated by adding some multiple of the standard deviation in
# some random (uniformly choosen) direction.
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
m2Angle=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  #p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,n)
  C[1:(p.out*n)] = 1
  C = sample(C)
  C = c(C, C)
  
  s=1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  angle = matrix(runif(n*p, 0, 2*pi),ncol=p)
  length = sapply(1:p, function(i){
    sample(size=n, sqrt(diag(omega))[(2*i-1):2*i], replace=T)
  } )
  length = matrix(length, ncol=p)
  k = t(rbind(cos(angle)*length, sin(angle)*length))*k
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away.  This function differs from m2Adj in that random
# errors (outliers) are generated by adding some multiple of the standard deviation in
# some random (uniformly choosen) direction.  However, the direction is restriced to be
# in the same direction as alpha, and hence outliers essentially increase the heaviness of
# the simulated tail.
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
#alphaWind: outliers must have errors added in the direction of alpha.  However, the angle
#  need not be alpha exactly, but instead alpha +/- alphaWind/2.  Defaults to 80 degrees.  
m2AngleRestrict=function(n,p,xi,omega,alpha,nu,p.out,k,alphaWind=2*pi/4.5){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  #p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,n)
  C[1:(p.out*n)] = 1
  C = sample(C)
  C = c(C, C)
  
  s=1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  length = NULL; angle = NULL
  for(i in 1:p){
    #Which values correspond to one pressure level?
    filt = (2*i-1):2*i
    length = cbind(length, sample(size=n, sqrt(diag(omega))[filt], replace=T) )
    alphaAngle = atan2(alpha[filt][1], alpha[filt][2])
    angle = cbind(angle, runif(n, alphaAngle-alphaWind/2, alphaAngle+alphaWind/2 ) )
  }
  length = matrix(length, ncol=p)
  angle = matrix(angle, ncol=p)
  k = t(rbind(cos(angle)*length, sin(angle)*length))*k
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#outlying at all higher levels
m3=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  
  C=rbinom(2*n,1,p.out)
  
  cout=which(C==1)
  nout=sum(C)
  s=2*rbinom(nout,1,0.5)-1
  cs.m=matrix(0,p,2*n,byrow=T)
  ti=runif(nout,1,p)
  
  for(j in 1:nout){
    part=ti[j]:p
    cs.m[part,cout[j]]=s[j]
  }
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}









#############################################################################################
#
#
# fastTLE implementations
#
#
#############################################################################################

#initial: vector for the initial guess for the parameters.  Results will depend on this initial
#   value, so it may be reasonable to run with several starting values.
#MLE: a function which takes two arguments: data and start.  data is supplied as an argument
#   to fastTLE.  This function should return a parameter vector of the same length as initial.
#   start may not do anything (as in the case of normal data) but must be an argument.  For
#   optimizations, start should be the best guess, and the last MLE will be supplied.
#negLogLik: a function which takes two arguments: param and data.  param should be a vector
#   of the parameters wished to optimize and data is supplied as an argument to fastTLE.  This
#   function is used to compute the negative log-likelihood of the observations and hence
#   determines which observations are used in computing the MLE at each iteration.  Thus, it
#   should return a vector of negative log-likelihoods for each observation.
#data: the data which should be fed into the data argument of likelihood.  data should be a
#   matrix of the data.  If it's not a matrix, fastTLE will attempt to coerce it to a matrix
#   (with a warning).
#k: A value between 0 and 1 specifying the % of observations that should be used to estimate
#   the MLE.  fastTLE will use the k*100 observations with the smallest likelihood.
#Note: This code was not designed to handle ties.  This will not be a problem with continuous
#   data, but there could be issues with discrete data.
fastTLE = function(initial, MLE, negLogLik, data, k, trace=F, ...){
  #Data quality checks
  if(!is(initial,"numeric"))
    stop("initial must be a numeric vector!")
  if(!is(MLE,"function"))
    stop("MLE must be a function!")
  if(!is(negLogLik,"function"))
    stop("MLE must be a function!")
  if(!"data" %in% names(formals(MLE)))
    stop("MLE must have a data argument!")
  if(any(names(formals(negLogLik))!=c("params","data")))
    stop("likelihood must have two arguments: params and data!")
  if(k>1 | k<=0)
    stop("k must be in (0,1]!")
  if(!is(data,"matrix")){
    warning("data is not a matrix!  Attempting to coerce...")
    data = matrix(data)
  }
  
  #Main implementation of code
  nUse = round(k*nrow(data))
  if(nUse<1)
    stop("nUse is smaller than 1.  k must be set to a larger value!")  
  if(trace) cat("Using", nUse, "observations for each MLE.\n")
  
  #Keep track of the obs used to fit the MLE each time.  Once they don't change, we've converged.
  swaps = 10 #Set it to anything >0 so while loop executes
  currMLE = initial
  if(trace) cat("Initial MLE:",initial,"\n")
  valsUsedOld = rep(F, nrow(data))
  iteration = 0
  #Algorithm will fail if negative log-likelihood increases.  So, set initial value to something huge.
  oldScores = rep(10000,nrow(data))
  while(swaps>0){
    scores = negLogLik(currMLE, data)
    if(sum(scores[valsUsed])>sum(oldScores[valsUsed]))
      stop("-log(likelihood(MLE; data)) increased.  Not possible with TLE.  Maybe MLE isn't converging?")
    else
      oldScores = scores
    
    valsUsed = rank(scores)<=nUse
    swaps = sum(valsUsed & !valsUsedOld)
    
    #Compute the MLE on the new dataset
    currMLE = MLE(data[valsUsed,,drop=F], start=currMLE)
    valsUsedOld = valsUsed
    iteration = iteration+1
    
    if(trace) cat("Iteration", iteration, "completed.  Current MLE:",round(currMLE,3),"Swaps:",swaps,"\n")
  }
  
  return(list(MLE=currMLE, valsUsed=valsUsed) )
}

#Wrapper to fastTLE assuming normal data.
fast.TLE.normal = function(data, k, initial=c(0,1)){
  #Data quality checks
  if(!is(data,"matrix")){
    warning("data is not a matrix.  Attempting to coerce...")
    data = matrix(data)
  }
  if(ncol(data)!=1)
    stop("This function is for univariate normal only!")
  if(k>1 | k<=0)
    stop("k must be in (0,1]")
  if(length(initial)!=2)
    stop("initial must be of length 2!")
  
  MLE = function(data, start){ c(mean(data), sd(data)*NROW(data)/(NROW(data)-1)) }
  negLogLik = function(params, data){
    mu = params[1]
    sigma = params[2]
    log(sigma) + (data-mu)^2/sigma^2
  }
  fastTLE(initial, MLE, negLogLik, data, k)
}

# data = rnorm(1000)
# data[1:100] = rnorm(100, sd=100)
# qplot(data)
# mean(data); sd(data)
# temp = fast.TLE.normal(data, k=.95)
# temp = fast.TLE.normal(data, k=.9)
# temp = fast.TLE.normal(data, k=.8)
# temp = fast.TLE.normal(data, k=.7)

#Wrapper to fastTLE assuming skew-t data.
fast.TLE.ST = function(data, k, initial=c(0,1,0,1000), ...){
  #Data quality checks
  if(!is(data,"matrix")){
    warning("data is not a matrix.  Attempting to coerce...")
    data = matrix(data)
  }
  if(ncol(data)!=1)
    stop("This function is for univariate skew-t only!")
  if(k>1 | k<=0)
    stop("k must be in (0,1]")
  if(length(initial)!=4)
    stop("initial must be of length 4!")
  
  MLE = function(data, start){ robustST(x=matrix(1,nrow(data)), y=data, robust=F, start=start) }
  negLogLik = function(params, data)
    sapply(data, st.pdev, dp=params, x=matrix(1))
  fastTLE(initial, MLE, negLogLik, data, k)
}

# data = rst(100, xi=0, omega=1, alpha=6, nu=12)
# qplot(data)
# data[1:10] = rnorm(10, sd=100)
# qplot(data)
# mean(data); sd(data)
# robustST(x=matrix(1,nrow=length(data)), y=data, robust=F)
# robustST(x=matrix(1,nrow=length(data)), y=data, robust=T)
# temp = fast.TLE.ST(data, k=.95)
# temp = fast.TLE.ST(data, k=.9)
# temp = fast.TLE.ST(data, k=.8)
# temp = fast.TLE.ST(data, k=.7)

#Wrapper to fastTLE assuming bivariate skew-t data.
fast.TLE.BST = function(data, k, initial=c(0,0,1,0,1,0,0,1000), ...){
  #Data quality checks
  if(!is(data,"matrix")){
    warning("data is not a matrix.  Attempting to coerce...")
    data = matrix(data)
  }
  if(ncol(data)!=2)
    stop("This function is for bivariate skew-t only!")
  if(k>1 | k<=0)
    stop("k must be in (0,1]")
  if(length(initial)!=8)
    stop("initial must be of length 8 (xi1, xi2, om11, om12, om22, al1, al2, nu)!")
  
  MLE = function(data, start){
    fit = robustST(x=matrix(1,nrow(data)), y=data, robust=F, start=start)
    return( c(fit$beta[1], fit$beta[2], fit$Omega[1,1], fit$Omega[1,2], fit$Omega[2,2]
              ,fit$alpha[1], fit$alpha[2], fit$nu) )
  }
  negLogLik = function(params, data){
    #We have to convert parameters to the type needed by mst.pdev:
    dp = list( beta=params[1:2], Omega=matrix(params[c(3,4,4,5)],nrow=2)
              ,alpha=params[6:7], nu=params[8])
    params = dplist2optpar(dp)

    sapply(1:nrow(data), function(i){
      mst.pdev(y=data[i,,drop=F], param=params, x=matrix(1))
    })
  }
  fit = fastTLE(initial, MLE, negLogLik, data, k)
  return( list(xi=fit$MLE[1:2], omega=matrix(fit$MLE[c(3,4,4,5)],nrow=2)
              ,alpha=fit$MLE[6:7], nu=fit$MLE[8]) )
}

# data = rmst(100, xi=c(0,0), Omega=diag(c(3,1)), alpha=c(6,-3), nu=12)
# qplot(data[,1], data[,2])
# data[1:3,] = rmst(3, Omega=diag(c(100,100)), alpha=c(6,-3))
# qplot(data[,1], data[,2])
# robustST(x=matrix(1,nrow=NROW(data)), y=data, robust=F)
# robustST(x=matrix(1,nrow=NROW(data)), y=data, robust=T)
# temp = fast.TLE.BST(data, k=.95)
# temp = fast.TLE.BST(data, k=.9)
# temp = fast.TLE.BST(data, k=.8)
# temp = fast.TLE.BST(data, k=.7)