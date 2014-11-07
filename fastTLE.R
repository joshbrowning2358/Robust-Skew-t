#initial: vector for the initial guess for the parameters.  Results will depend on this initial
#   value, so it may be reasonable to run with several starting values.
#MLE: a function which takes one argument: data.  data is supplied as an argument to fastTLE.
#   This function should return a parameter vector of the same length as initial.  May
#   accept additional arguments specified via ...
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
fastTLE = function(initial, MLE, negLogLik, data, k, ...){
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
  cat("Using", nUse, "observations for each MLE.\n")
  
  #Keep track of the obs used to fit the MLE each time.  Once they don't change, we've converged.
  swaps = 10 #Set it to anything >0 so while loop executes
  currMLE = initial
  cat("Initial MLE:",initial,"\n")
  valsUsedOld = rep(F, nrow(data))
  iteration = 0
  while(swaps>0){
    scores = negLogLik(currMLE, data)
    valsUsed = rank(scores)<=nUse
    swaps = sum(valsUsed & !valsUsedOld)
    
    #Compute the MLE on the new dataset
    currMLE = MLE(data[valsUsed,])
    valsUsedOld = valsUsed
    iteration = iteration+1
    
    cat("Iteration", iteration, "completed.  Current MLE:",round(currMLE,3),"Swaps:",swaps,"\n")
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
  
  MLE = function(data){ c(mean(data), sd(data)*NROW(data)/(NROW(data)-1)) }
  negLogLik = function(params, data){
    mu = params[1]
    sigma = params[2]
    log(sigma) + (data-mu)^2/sigma^2
  }
  fastTLE(initial, MLE, negLogLik, data, k)
}

data = rnorm(1000)
data[1:100] = rnorm(100, sd=100)
qplot(data)
mean(data); sd(data)
temp = fast.TLE.normal(data, k=.95)
temp = fast.TLE.normal(data, k=.9)
temp = fast.TLE.normal(data, k=.8)
temp = fast.TLE.normal(data, k=.7)

#Wrapper to fastTLE assuming normal data.
fast.TLE.skewt = function(data, k, initial=c(0,1,0,1000)){
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
  
  MLE = function(data){ c(mean(data), sd(data)*NROW(data)/(NROW(data)-1)) }
  negLogLik = function(params, data){
    mu = params[1]
    sigma = params[2]
    log(sigma) + (data-mu)^2/sigma^2
  }
  fastTLE(initial, MLE, negLogLik, data, k)
}