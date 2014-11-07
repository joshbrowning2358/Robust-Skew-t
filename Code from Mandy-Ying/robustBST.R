#setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
library(sn)
source("sn-funct.R")
source("modified_mst.pdev.grad.R")

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
    return( list(beta=optpar$beta, Omega=optpar$Omega, alpha=optpar$alpha, df=optpar$nu, convergence=fit[length(fit)] ) )
  }
}
