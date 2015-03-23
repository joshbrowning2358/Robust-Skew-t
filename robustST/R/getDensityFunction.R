##' Get Density Functions
##' 
##' In order to fit different robust distributions, we need to use different
##' density and gradient of density functions.  This function contains the
##' function definitions for four classes of families:  N (normal), T
##' (Student's-t), SN (skew normal), ST (skew t).
##' 
##' @param family The family to use, should be "N", "T", "SN" or "ST".
##' @param dimension How many columns does y have?
##' @param robust A logical value indicating whether the density and gradient
##' returned should correspond to the robust estimators or not.
##' 
##' @return A list of two functions, the first being the density function and
##' the second being the gradient of the density function.  To understand the
##' usage of this function, see ?robustSTOnceK.
##' 

## Test that this function works for all distributions

getDensityFunction = function(family, dimension, robust = TRUE){
    ## Data quality checks
    if(! family %in% c("N", "T", "SN", "ST"))
        stop("Invalid family provided.  Must be 'N', 'T', 'SN' or 'ST'")
    
    if(family == "N"){
        if(robust){
            density = function(param, y, k){
                p = NCOL(y)
                ## Add alpha and nu to param vector
                param = c(param, rep(0, p), Inf)
                mst.pdev.robust(param, y, k, symmetr = TRUE, fixed.nu = TRUE)
            }
            gradient = function(param, y, k){
                p = NCOL(y)
                ## Add alpha and nu to param vector
                param = c(param, rep(0, p), Inf)
                mst.pdev.grad.robust(param, y, k, symmetr = TRUE,
                                     fixed.nu = TRUE)
            }
        } else {
            density = function(param, y, k){
                p = NCOL(y)
                ## Add alpha and nu to param vector
                param = c(param, rep(0, p), Inf)
                sn:::mst.pdev(param, x = matrix(1, nrow = NROW(y)), y, w,
                              symmetr = TRUE, fixed.nu = TRUE)
            }
            gradient = function(param, y, w = rep(1, NROW(y))){
                p = NCOL(y)
                ## Add alpha and nu to param vector
                param = c(param, rep(0, p), Inf)
                sn:::mst.pdev.grad(param, x = matrix(1, nrow = NROW(y)), y, w,
                                   symmetr = TRUE, fixed.nu = TRUE)
            }
        }
    } else if(family == "T"){
        if(robust){
            density = function(param, y, k){
                p = NCOL(y)
                ## Add alpha to param vector
                param = c(param[-length(param)], rep(0, p),
                          param[length(param)])
                mst.pdev.robust(param, y, k, symmetr = TRUE)
            }
            gradient = function(param, y, k){
                p = NCOL(y)
                ## Add alpha to param vector
                param = c(param[-length(param)], rep(0, p),
                          param[length(param)])
                mst.pdev.grad.robust(param, y, k, symmetr = TRUE)
            }
        } else {
            density = function(param, y, k){
                p = NCOL(y)
                ## Add alpha to param vector
                param = c(param[-length(param)], rep(0, p),
                          param[length(param)])
                sn:::mst.pdev(param, x = matrix(1, nrow = NROW(y)), y, w,
                              symmetr = TRUE)
            }
            gradient = function(param, y, w = rep(1, NROW(y))){
                p = NCOL(y)
                ## Add alpha to param vector
                param = c(param[-length(param)], rep(0, p),
                          param[length(param)])
                sn:::mst.pdev.grad(param, x = matrix(1, nrow = NROW(y)), y, w,
                                   symmetr = TRUE)
            }
        }
    } else if(family == "SN"){
        if(robust){
            density = function(param, y, k){
                param = c(param, Inf)
                mst.pdev.robust(param, y, k, fixed.nu = Inf)
            }
            gradient = function(param, y, k){
                param = c(param, Inf)
                mst.pdev.grad.robust(param, y, k, fixed.nu = Inf)
            }
        } else {
            density = function(param, y, k){
                param = c(param, Inf)
                sn:::mst.pdev(param, x = matrix(1, nrow = NROW(y)), y, w,
                              fixed.nu = Inf)
            }
            gradient = function(param, y, w = rep(1, NROW(y))){
                param = c(param, Inf)
                sn:::mst.pdev.grad(param, x = matrix(1, nrow = NROW(y)), y, w,
                                   fixed.nu = Inf)
            }
        }
    } else if(family == "ST"){
        if(robust)
            return(list(density = mst.pdev.robust,
                        gradient = mst.pdev.grad.robust))
        else
            return(list(density = sn:::mst.pdev,
                        gradient = sn:::mst.pdev.grad))
    } else {
        stop("Current family not yet implemented")
    }
    stop("Something weird happened.  This line should never be reached.")
}