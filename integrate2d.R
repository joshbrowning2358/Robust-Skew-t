library(sn)
library(reshape)

xi = runif(2, -2, 2)
omega = matrix(c(1,0,0,1), nrow=2)
alpha = runif(2, -5, 5)
nu = exp(runif(1,log(4),8))
param = list(xi, omega, alpha, nu)

xi2 = xi + rnorm(2)
omega2 = omega
alpha2 = alpha + rnorm(2)
nu2 = nu*exp(rnorm(1))
param2 = list(xi2, omega2, alpha2, nu2)

diffEst = function(param, param2, xGridPts=30, yGridPts=xGridPts){
  grid = defineGrid(param, xGridPts, yGridPts)
  int = integrate2d(grid[[1]], grid[[2]], function(x){
    abs(dmst(x, xi=param[[1]], Omega=param[[2]], alpha=param[[3]], nu=param[[4]]) -
      dmst(x, xi=param2[[1]], Omega=param2[[2]], alpha=param2[[3]], nu=param2[[4]]))
  })
  return(int)
}

defineGrid = function(param, xGridPts=30, yGridPts=xGridPts){
  xi = param[[1]]
  omega = param[[2]]
  alpha = param[[3]]
  nu = param[[4]]

  #Naive approach: assume bivariate skew-t is best represented by two marginal skew-t's.
  #xGrid = qst( 1:(xGridPts-1)/xGridPts, xi[1], omega[1,1], alpha[1], nu )
  #yGrid = qst( 1:(yGridPts-1)/yGridPts, xi[2], omega[2,2], alpha[2], nu )
  #return(list(xGrid, yGrid))
  
  #Naive approach 2: Step by omega and go xGridPts/2*omega or yGridPts/2*omega in +/- directions
  xGrid = param[[1]][1] + (-xGridPts/2):(xGridPts/2)*param[[2]][1,1]/2
  yGrid = param[[1]][2] + (-yGridPts/2):(yGridPts/2)*param[[2]][2,2]/2
  return(list(xGrid, yGrid))
}

#Integrates func() over x and y space.  The integral is estimated using a 2-d analogue of the
# trapesoidal rule.  Each "trapesoid" has a rectangular base and 4 (possibly different)
# heights of the corners.  We compute the volume as A(z1+z2+z3+z4)/4, where A is the base area
# and zi the ith height.  Note that trapesoids are only defined within the grid, thus all
# volumne outside of the grid is ignored.
#xGrid: the x values to approximate the integral over.
#yGrid: the y values to approximate the integral over.
#func: a function whose first argument is a 2d vector, and the function evaluates to a scalar.
#...: additional arguments to pass to func.
integrate2d = function( xGrid, yGrid, func, ...){
  library(reshape)
  
  grid = merge( xGrid, yGrid)
  colnames(grid) = c("x", "y")

  vals = apply( grid, 1, func, ... )
  data = melt( cbind(grid, vals), id.vars=1:3 )
  data = cast( data, x ~ y, value="vals")
  #Remove xGrid from data (just want data values in the matrix)
  rownames(data) = data[,1]
  data[,1] = NULL
  
  #Average columns of data
  data2 = matrix(0, nrow=nrow(data), ncol=nrow(data)-1)
  for( i in 1:(ncol(data)-1) ){
    data2[,i] = (data[,i] + data[,i+1])/2
  }
  
  data = data2
  data2 = matrix(0, nrow(data)-1, ncol=ncol(data) )
  for( i in 1:(nrow(data)-1) ){
    data2[i,] = (data[i,] + data[i+1,])/2
  }
  
  widths = diff( xGrid )
  heights = diff( yGrid )
  areas = t(t(widths))%*%t(heights)
  volumes = areas*data2
  return(sum(volumes))
}

#Test it out
#integrate2d( xGrid=-2:2, yGrid=-2:2, func=function(x){dmst(x, Omega=diag(c(1,1)), alpha=c(0,0))} )
#integrate2d( xGrid=-5:5, yGrid=-5:5, func=function(x){dmst(x, Omega=diag(c(1,1)), alpha=c(0,0))} )
n = 10
#integrate2d( xGrid=qnorm(1:(n-1)/n), yGrid=qnorm(1:(n-1)/n), func=function(x){dmst(x, Omega=diag(c(1,1)), alpha=c(0,0))} )
integrate2d( xGrid=(-n/2):(n/2), yGrid=(-n/2):(n/2), func=function(x){dmst(x, xi=c(0,0), Omega=diag(c(1,1)), alpha=c(0,0))} )
