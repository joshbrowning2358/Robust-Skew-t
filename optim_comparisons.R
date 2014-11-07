library(sn)
library(rbenchmark)

set.seed(123)
d = rnorm(10000)
nLL = function(pars){
  xi = pars[1]
  omega = pars[2]
  alpha = pars[3]
  nu = pars[4]
  sum( -log( dst(d, xi, omega, alpha, nu) ) )
}

#Does not converge at 20 iterations
optim( c(4, 2, 1, 100), nLL, lower=c(-Inf,0,-Inf,0), method="L-BFGS-B", control=list(maxit=20, trace=T) )

benchmark(
  Once = optim( c(4, 2, 1, 100), nLL, lower=c(-Inf,0,-Inf,0), method="L-BFGS-B"
      ,control=list(maxit=20) )
  ,Iter = {
    par = c(4,2,1,100)
    for(i in 1:20){
    par = optim( par, nLL, lower=c(-Inf,0,-Inf,0), method="L-BFGS-B"
      ,control=list(maxit=1) )$par
    }
  }
  ,replications=3 )
