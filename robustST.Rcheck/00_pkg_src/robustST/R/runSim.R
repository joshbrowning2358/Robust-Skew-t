##' Run simulation
##' 
##' This function runs a simulation to measure the performance of the robust
##' and non-robust fitting methods on data where the distribution is known.
##' 
##' @param n Sample size to simulate.
##' @param outPct Percent of outliers.  Actual number is round(n*outPct)
##' @param outSigma: Outliers are created by adding on an error of
##' N(0,outSigma^2)
##' @param k A numeric vector of constants to use in "capping" the likelihood
##' function.  Each constant is tried in turn, and results are returned for all
##' values.
##' @param fast.k A vector of k values for the fastTLE algorithm.  Each value
##' should represent the proportion of data to be used in estimating the MLE.
##' @param xi0 Simulated center parameter of the skew-t to simulate.
##' @param omega0 Simulated scale parameter of the skew-t to simulate.
##' @param alpha0 Simulated skewness parameter of the skew-t to simulate.
##' @param nu0 Simulated heaviness of tails parameter of the skew-t to simulate.
##' @param pressure Use Denver station data at this pressure level, and
##' simulate the skew-t using Ying's code
##' @param type The type of skewness to use.  Must be "MVN", "obs", or "EX".
##' @param restrict Should outliers be generated in the tail of the skew-t
##' only?  Only applies if pressure and type are not null.
##' 
##' @return A data.frame containing the results of the simulations.
##' 
##' @export
##' 

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
            y = as.vector(sn:::rst(n, xi0, Omega0, alpha0, nu0))    
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
