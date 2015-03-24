##' Robust adjustment of the negative log-likelihood
##' 
##' Function implementing the robustness adjustment to the negative
##' log-likelihood.
##' 
##' Essentially, if the negative log-likelihood value is too large, it is
##' bounded via Psi().  The psi() function bounds the derivative.
##' 
##' @param originalNLL The original negative log-likelihood
##' @param k The parameter controlling the amount of robustification.  Negative
##' log-likelihoods larger than k are reduced, and the adjusted value is always
##' less than k+1.
##' 
##' @return The adjusted negative log-likelihood.
##' 
##' @export
##' 

Psi = function(originalNLL, k){
    ifelse(originalNLL > k,
           k + (1 - exp(-(originalNLL - k))),
           originalNLL)
}