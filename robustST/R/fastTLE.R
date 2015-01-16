##' Fast TLE
##' 
##' An method of generating robust estimates for parameters of a distribution,
##' different from the capped likelihood approach.  Rather than capping the
##' deviance, a few of the worst observations (i.e. largest deviance) are
##' not used in the estimation of the MLE.  Thus, this procedure is iterative;
##' the worst observations may change as the MLE changes, and the MLE may
##' change as the set of observations used changes.  Fast TLE is an algorithm
##' which is NOT guaranteed to find the optimum, but typically obtains a good
##' approximation quickly.
##' 
##' @param initial A vector for the initial guess for the parameters.  Results
##' will depend on this initial value, so it may be reasonable to run with
##' several starting values.
##' @param MLE A function which takes two arguments: data and start.  data is
##' supplied as an argument to fastTLE.  This function should return a
##' parameter vector of the same length as initial. start may not do anything
##' (as in the case of normal data) but must be an argument.  For
##' optimizations, start should be the best guess, and the last MLE will be
##' supplied.
##' @param negLogLik A function which takes two arguments: param and data.
##' param should be a vector of the parameters wished to optimize and data is
##' supplied as an argument to fastTLE.  This function is used to compute the
##' negative log-likelihood of the observations and hence determines which
##' observations are used in computing the MLE at each iteration.  Thus, it
##' should return a vector of negative log-likelihoods for each observation.
##' @param data The data which should be fed into the data argument of
##' likelihood.  data should be a matrix of the data.  If it's not a matrix,
##' fastTLE will attempt to coerce it to a matrix (with a warning).
##' @param k A value between 0 and 1 specifying the % of observations that
##' should be used to estimate the MLE.  fastTLE will use the k*100
##' observations with the smallest likelihood.  Note: This code was not
##' designed to handle ties.  This will not be a problem with continuous data,
##' but there could be issues with discrete data.
##' @param trace Should output be printed as the algorithm proceeds?
##' @param ... Currently not implemented.
##' 
##' @return A list of two elements: the estimated MLE and a vector specifying
##' which observations were used in the estimation.
##' 
##' @export
##' 

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