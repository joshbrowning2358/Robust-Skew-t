mst.pdev.grad.new <- function(param, x, y, w, fixed.nu=NULL, penalty=NULL, trace=FALSE)
{
  #Initialization
  d <- ncol(y)
  p <- ncol(x)

  #Extract parameters from param
  beta<- matrix(param[1:(p*d)],p,d)
  D  <- exp(-2*param[(p*d+1):(p*d+d)])
  A  <- diag(d)
  i0 <- p*d+d*(d+1)/2
  if(d>1) A[!lower.tri(A,diag=TRUE)] <- param[(p*d+d+1):i0]
  eta   <- param[(i0+1):(i0+d)]
  nu <- if(is.null(fixed.nu))  exp(param[i0+d+1]) else fixed.nu

  #Define additional parameters
  Oinv  <- t(A) %*% diag(D,d,d) %*% A
  u     <- y - x %*% beta
  Q     <- as.vector(apply((u %*% Oinv)*u,1,sum))
  L     <- as.vector(u %*% eta)
  sf    <- if(nu<10000) sqrt((nu+d)/(Q+nu)) else sqrt((1+d/nu)/(1+Q/nu))
  t.    <- L*sf

  dlogft<- (-0.5)*(1+d/nu)/(1+Q/nu)
  dt.dL <- sf
  dt.dQ <- (-0.5)*L*sf/(Q+nu)
  logT. <- pt(t., nu+d, log.p=TRUE)
  dlogT.<- exp(dt(t., nu+d, log=TRUE) - logT.)
  u.w<- u*w
  Dbeta <- (-2* t(x) %*% (u.w*dlogft) %*% Oinv 
            - outer(as.vector(t(x) %*% (dlogT. * dt.dL* w)), eta)
            - 2* t(x) %*% (dlogT.* dt.dQ * u.w) %*% Oinv )
  Deta  <- apply(dlogT.*sf*u.w, 2, sum)
  if(d>1) {
     M  <- 2*( diag(D,d,d) %*% A %*% t(u * dlogft
               + u * dlogT. * dt.dQ) %*% u.w)
     DA <- M[!lower.tri(M,diag=TRUE)]
     }
  else DA<- NULL
#  M     <- ( A %*% t(u*dlogft + u*dlogT.*dt.dQ) %*% u.w %*% t(A))
#  if(d>1) DD <- diag(M) + 0.5*sum(w)/D
#     else DD <- as.vector(M + 0.5*sum(w)/D) 
#  grad <- (-2)*c(Dbeta,DD*(-2*D),DA,Deta)

  myDD = -NROW(y) + t(t(dlogft)%*%t(A%*%t(u))^2)*(-2*D) + t((dlogT.*dt.dQ)%*%t(A%*%t(u))^2)*(-2*D)
  grad <- (-2)*c(Dbeta,myDD,DA,Deta)

  if(is.null(fixed.nu)) {
    df0 <- min(nu, 1e8)
    if(df0 < 10000){
       diff.digamma <- digamma((df0+d)/2) - digamma(df0/2)
       log1Q<- log(1+Q/df0)
     }
    else
      {
       diff.digamma <- log1p(d/df0)
       log1Q <- log1p(Q/df0)
      }
    dlogft.ddf <- 0.5 * (diff.digamma - d/df0
                        + (1+d/df0)*Q/((1+Q/df0)*df0) - log1Q)
    eps   <- 1.0e-4
    df1 <- df0 + eps
    sf1 <- if(df0 < 1e4) sqrt((df1+d)/(Q+df1)) else sqrt((1+d/df1)/(1+Q/df1))
    logT.eps <- pt(L*sf1, df1+d, log.p=TRUE)
    dlogT.ddf <- (logT.eps-logT.)/eps
    Ddf   <- sum((dlogft.ddf + dlogT.ddf)*w)
    grad <- c(grad, -2*Ddf*df0)
    }
  if(!is.null(penalty)) { 
    Ainv <- backsolve(A, diag(d))
    Omega <- Ainv %*% diag(1/D,d,d) %*% t(Ainv)
    omega <- diag(Omega)
    alpha <- eta*omega
    Q <- Qpenalty(list(alpha, cov2cor(Omega)), nu, der=1)
    comp <-  1:(length(alpha)+is.null(fixed.nu))
    Qder <- attr(Q, "der1") * c(1/omega, 1)[comp] 
    # gradient for transformed variable (alpha --> eta)
    grad <- grad + 2*c(rep(0, p*d + d*(d+1)/2),  Qder)
    }
  if(trace) cat("mst.pdev.grad: norm is ", format(sqrt(sum(grad^2))), "\n")  
  return(grad)
}