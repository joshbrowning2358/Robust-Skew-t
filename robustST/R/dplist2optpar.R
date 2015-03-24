dplist2optpar <- function (dp, Omega.inv = NULL) 
{
    beta <- dp[[1]]
    Omega <- as.matrix(dp[[2]])
    alpha <- dp[[3]]
    d <- length(alpha)
    nu <- if (is.null(dp$nu)) 
        NULL
    else dp$nu
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
    return(param)
}