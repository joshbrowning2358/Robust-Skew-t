##' Get Dimension from Parameter Vector Size
##' 
##' The parameter vector size is directly related to the dimension of the
##' dependent data, but this relationship isn't directly obvious.  For example,
##' for 1-d, we have 4 parameters (xi, omega, alpha, nu).  For 2-d, we have 8
##' parameters (2 for xi and alpha, 3 for Omega, and 1 for nu).  For 3-d, we
##' have 13 parameters (3 for xi and alpha, 6 for Omega, and 1 for nu).  This
##' relationship can be determined from the quadratic formula
##' (d + d(d+1)/2 + d + 1 = length)
##' and the simple solution to this equation is provided in this function.
##' 
##' @param length The length of the parameter vector, expressed as a single
##' numerical value.
##' 
##' @return The dimension of the data (a single integer value).
##' 

dimensionFromParamSize = function(length){
    out = -2.5 + .5 * sqrt(17 + 8 * length)
    if(out != floor(out))
        stop("The length of the parameter vector isn't right, since the ",
             "answer is not an integer.")
    return(out)
}
