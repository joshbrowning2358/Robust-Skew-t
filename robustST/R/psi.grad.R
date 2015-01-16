##' Robust adjustment of the derivative of the negative log-likelihood
##' 
##' This function implements an adjustment to the derivative of the density
##' function to make estimation of parameters more robust.  If a negative
##' log-likelihood is supplied, then psi returns an adjusted negative
##' log-likelihood. Essentially, if the negative log-likelihood value is too
##' large, it is bounded via this function.
##' 
##' Usually, we're interested in df/dp, where f is the negative log-likelihood
##' and p a parameter.  df/dp should be updated with df/dp*psi(f, p) to be
##' "robustified".
##' 
##' @param originalNLL The original negative log-likelihood
##' @param k A parameter controlling the amount of robustification.  Negative
##' log-likelihoods larger than k are reduced, and the adjusted value is always
##' less than 2*k.
##' 
##' @return Returns 
##' 
##' @export
##' 

psi.grad = function(originalNLL, k){
    ifelse(originalNLL>k,
           exp(-originalNLL/k+1),
           1)
}
