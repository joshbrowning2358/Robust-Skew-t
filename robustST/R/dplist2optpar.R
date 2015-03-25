dplist2optpar <- function (dp, family, Omega.inv = NULL) 
{
    
    ## Data Quality Checks
    stopifnot(family %in% c("N", "T", "SN", "ST"))
    stopifnot(names(dp) %in% c("xi", "Omega", "alpha", "nu"))
    
    beta <- dp$xi
    Omega <- as.matrix(dp$Omega)
    if(family %in% c("N", "T"))
        dp$alpha = rep(0, length(beta))
    if(family %in% c("N", "SN"))
        dp$nu = Inf
    alpha <- dp$alpha
    d <- length(alpha)
    nu = dp$nu
    eta <- alpha/sqrt(diag(Omega))
    Oinv <- if (is.null(Omega.inv)) 
        mnormt:::pd.solve(Omega)
    else Omega.inv
    if (is.null(Oinv)) 
        stop("matrix Omega not symmetric positive definite")
    upper <- chol(Oinv)
    D <- diag(upper)
    A <- upper/D
    D <- D^2
    param <- if (d > 1) 
        c(beta, -log(D)/2, A[!lower.tri(A, diag = TRUE)], eta)
    else c(beta, -log(D)/2, eta)
    if (!is.null(dp$nu)) 
        param <- c(param, log(dp$nu))
    param <- as.numeric(param)
    attr(param, "ind") <- cumsum(c(length(beta), d, d * (d - 
        1)/2, d, length(dp$nu)))
    if(family == "N")
        param = param[1:(d+d*(d+1)/2)]
    if(family == "T")
        param = c(param[1:(d+d*(d+1)/2)], param[length(param)])
    if(family == "SN")
        param = param[-length(param)]
    return(param)
}