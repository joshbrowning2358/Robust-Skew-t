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

getDensityFunction = function(family, dimension, robust = TRUE){
    ## Data quality checks
    if(! family %in% c("N", "T", "SN", "ST"))
        stop("Invalid family provided.  Must be 'N', 'T', 'SN' or 'ST'")
    
    if(family == "N"){
        if(robust & dimension == 1)
            return(list(density = n.dev.robust,
                        gradient = n.dev.gh.robust))
        else if(!robust & dimension == 1)
            return(list(density = n.dev,
                        gradient = n.dev.gh))
        stop("Current family/inputs not yet implemented")                            
    } else if(family == "T"){
        stop("Current family not yet implemented")                            
    } else if(family == "SN"){
        if(!robust & dimension == 1)
            return(list(density = sn:::sn.pdev,
                        gradient = sn:::sn.pdev.gh))
        else if(!robust & dimension > 1)
            return(list(density = sn:::msn.dev,
                        gradient = sn:::msn.dev.grad))
        stop("Current family not yet implemented")                            
    } else if(family == "ST"){
        if(dimension == 1 & robust)
            return(list(density = st.pdev.robust,
                        gradient = st.pdev.gh.robust))
        else if(dimension > 1 & robust)
            return(list(density = mst.pdev.robust,
                        gradient = mst.pdev.grad.robust))
        else if(dimension == 1 & !robust)
            return(list(density = sn:::st.pdev,
                        gradient = sn:::st.pdev.grad))
        else if(dimension > 1 & !robust)
            return(list(density = sn:::st.pdev,
                        gradient = sn:::st.pdev.gh))
    } else {
        stop("Current family not yet implemented")
    }
    stop("Something weird happened.  This line should never be reached.")
}