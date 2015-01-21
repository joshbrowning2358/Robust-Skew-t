#setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
source("Code/sn-funct.R")
source("Code/modified_mst.pdev.grad.R")

#Original mst.pdev.grad is not additive for D matrix.  New implementation fixes this.
mst.pdev.grad = mst.pdev.grad.new





#This function approximates the integral of the abs. value of the difference in 2d skew-t densities.
diffEst = function(param, param2, xGridPts=30, yGridPts=xGridPts){
  grid = defineGrid(param, xGridPts, yGridPts)
  int = integrate2d(grid[[1]], grid[[2]], function(x){
    abs(dmst(x, xi=param[[1]], Omega=param[[2]], alpha=param[[3]], nu=param[[4]]) -
      dmst(x, xi=param2[[1]], Omega=param2[[2]], alpha=param2[[3]], nu=param2[[4]]))
  })
  return(int)
}

#This function chooses a grid for evaluating the 2d integral of the difference in densities.
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

#From ying_model3.R (Denver station parameters)
denverParams = function(type=c("MVN", "obs", "EX")){
  if(length(type)>1)
    type==type[1]
  if(!type %in% c("MVN", "obs", "EX"))
    stop("type must be one of 'MVN', 'obs', or 'EX'!")
  
  xi=c(3.5732607, 10.3322993, 13.0043508, 16.5135646, 19.0336571, 20.8985299, 11.4608094, 5.2606013, -1.0291708, -2.2748794, -1.7431219, -1.2281764, -0.9888074, -0.4555666, -0.7184110, -0.9835646)
  names(xi)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  
  omega=rbind( c(22.3972699, 13.984818,  16.039642,  17.684241,  18.63175,  18.346261, 10.4795136,  8.2677082, -3.9072388, -10.004767, -10.995055, -13.833277, -13.634031, -10.2153270, -2.815235, -0.5339653), c(13.9848178, 53.727990,  59.020103,  66.056794,  66.75493,  59.537442, 23.6810410, 14.8809412,  5.5295733,   3.477148,   6.406637,   8.474633,   9.738715,  11.3172025,  4.378073,  3.3395517), c(16.0396421, 59.020103,  90.528398, 102.709659, 103.20949,  88.314971, 30.7581432, 17.6257511,  9.3983732,  11.193745,  17.330145,  24.295520,  26.207537,  25.1415773,  8.529959,  6.3410444), c(17.6842407, 66.056794, 102.709659, 147.899721, 148.89594, 124.512607, 38.0446733, 19.9519976, 12.4206720,  18.832824,  27.023667,  38.891596,  42.174570,  38.6743609, 12.975820,  8.9548713), c(18.6317507, 66.754935, 103.209486, 148.895942, 171.79385, 145.080702, 43.3933979, 22.0252564, 12.3625644,  19.464066,  27.533566,  40.567982,  44.611103,  40.2599930, 13.941600,  9.4114303), c(18.3462610, 59.537442,  88.314971, 124.512607, 145.08070, 149.877677, 46.6103745, 24.3935571, 10.9785169,  14.172450,  21.046724,  30.593732,  35.229002,  33.2482152, 10.321343,  7.7821568), c(10.4795136, 23.681041,  30.758143,  38.044673,  43.39340,  46.610374, 51.6696619, 33.6064021,  0.1524421,  -4.419451,  -4.078705,  -3.998690,  -3.327270,  -0.0476168, -1.178374,  1.0679616), c( 8.2677082, 14.880941,  17.625751,  19.951998,  22.02526,  24.393557, 33.6064021, 40.3398453, -2.7411404,  -7.184680,  -7.805223,  -9.335120,  -9.578368,  -5.4294290, -2.678479, -0.6314295), c(-3.9072388,  5.529573,   9.398373,  12.420672,  12.36256,  10.978517,  0.1524421, -2.7411404, 25.9163894,  27.043738,  30.793393,  37.083542,  38.886431,  35.2738880, 13.268837,  6.7957864),c(-10.0047669,  3.477148,  11.193745,  18.832824,  19.46407,  14.172450, -4.4194513, -7.1846796, 27.0437381,  78.522421,  90.324931, 108.070579, 112.710362, 100.0875035, 36.757158, 17.5988472),c(-10.9950547,  6.406637,  17.330145,  27.023667,  27.53357,  21.046724, -4.0787046, -7.8052225, 30.7933930,  90.324931, 126.802112, 151.111722, 155.933320, 135.8397389, 47.631396, 22.8484928),c(-13.8332773,  8.474633,  24.295520,  38.891596,  40.56798,  30.593732, -3.9986897, -9.3351204, 37.0835417, 108.070579, 151.111722, 209.976779, 215.636087, 185.3566163, 62.552556, 29.9704567),c(-13.6340314,  9.738715,  26.207537,  42.174570,  44.61110,  35.229002, -3.3272702, -9.5783683, 38.8864311, 112.710362, 155.933320, 215.636087, 242.728801, 208.9461148, 70.680936, 33.6640122),c(-10.2153270, 11.317202,  25.141577,  38.674361,  40.25999,  33.248215, -0.0476168, -5.4294290, 35.2738880, 100.087504, 135.839739, 185.356616, 208.946115, 208.5533958, 72.232734, 34.7282020), c(-2.8152352,  4.378073,   8.529959,  12.975820,  13.94160,  10.321343, -1.1783744, -2.6784786, 13.2688371,  36.757158,  47.631396,  62.552556,  70.680936,  72.2327344, 46.452689, 21.1557030), c(-0.5339653,  3.339552,   6.341044,   8.954871,   9.41143,   7.782157,  1.0679616, -0.6314295,  6.7957864,  17.598847,  22.848493,  29.970457,  33.664012,  34.7282020, 21.155703, 18.8168900))
  
  colnames(omega)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  rownames(omega)=c("u700", "u500", "u400", "u300", "u250", "u200", "u100", "u70",  "v700", "v500", "v400", "v300", "v250", "v200", "v100", "v70" )
  
  #Constructing an omega that is perfectly symmetric.  omega above fails the symmetry test because of round-off error.
  UT=upper.tri(omega)
  omega2=matrix(0,nrow=nrow(omega),ncol=ncol(omega))
  omega2[UT]=omega[UT]
  omega2=omega2+t(omega2)
  diag(omega2)=diag(omega)
  eigen(omega2)$values
  colnames(omega2) = colnames(omega)
  rownames(omega2) = rownames(omega)
  
  #Values based on those observed at Denver Station
  alpha.obs=c(2.16,  1.44,  1.35,  0.95,  0.60,  0.31,  2.62,  3.02, -0.03, -0.81, -0.82, -0.89, -0.82, -0.75, -0.31,  0.39 )
  #Rewrite so we have a 1,-1 for the 300 pressure levels
  alpha.obs=c(2.16,  1.44,  1.35,  1,  0.60,  0.31,  2.62,  3.02, -0.03, -0.81, -0.82, -1, -0.82, -0.75, -0.31,  0.39 )
  #alpha.obs=c(0.09, -0.07,  0.13, -0.25,  0.18, -0.08, -0.04,  0.00, -0.09,  0.11,  0.02, -0.08, -0.03,  0.07,  0.07, -0.10)
  df.obs=10
  names(alpha.obs) = names(xi)

  #Values for a MVN distribution
  alpha.MVN=rep(0,16)
  df.MVN=Inf
  names(alpha.MVN) = names(xi)
  
  #Values for distributions that are more skewed than what was observed at the Denver Station
  alpha.EX=c(3:10,3:10)
  #Rewrite so we have a 6,0 for the 300 pressure levels
  alpha.EX[4]=6; alpha.EX[12]=0
  df.EX=5
  names(alpha.EX) = names(xi)
  
  if(type=="MVN")
    return(list( xi=xi, omega=omega2, alpha=alpha.MVN, nu=df.MVN) )
  if(type=="obs")
    return(list( xi=xi, omega=omega2, alpha=alpha.obs, nu=df.obs) )
  if(type=="EX")
    return(list( xi=xi, omega=omega2, alpha=alpha.EX, nu=df.EX) )
}








#############################################################################################
#
#
# Data simulation functions from Ying
#
#
#############################################################################################

#model 1 without outliers
m1=function(n,p,xi,omega,alpha,nu){
  #p by n data matrix
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  list(u=u,v=v)
}

#outlying at all levels
m2=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  C=rbinom(2*n,1,p.out)
  
  s=2*rbinom(2*n,1,0.5)-1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
m2Adj=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,2*n) #JB
  C[1:(p.out*n)] = 1 #JB
  C = sample(C) #JB
  
  s=2*rbinom(2*n,1,0.5)-1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  k = sqrt(rep( diag(omega), each=n))*k #JB
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away.  This function differs from m2Adj in that random
# errors (outliers) are generated by adding some multiple of the standard deviation in
# some random (uniformly choosen) direction.
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
m2Angle=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  #p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,n)
  C[1:(p.out*n)] = 1
  C = sample(C)
  C = c(C, C)
  
  s=1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  angle = matrix(runif(n*p, 0, 2*pi),ncol=p)
  length = sapply(1:p, function(i){
    sample(size=n, sqrt(diag(omega))[(2*i-1):2*i], replace=T)
  } )
  length = matrix(length, ncol=p)
  k = t(rbind(cos(angle)*length, sin(angle)*length))*k
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#Have k control the number of sd away.  This function differs from m2Adj in that random
# errors (outliers) are generated by adding some multiple of the standard deviation in
# some random (uniformly choosen) direction.  However, the direction is restriced to be
# in the same direction as alpha, and hence outliers essentially increase the heaviness of
# the simulated tail.
#n: the number of observations to simulate
#p: the number of pressure levels simulated.  length(xi)=length(alpha)=dim(omega)=p*2
#xi, omega, alpha, nu: skew-t parameters.  Typically estimated params from Denver station.
#p.out: Originally probability of outlier.  Updated to proportion of outliers.
#k: Originally constant to add/subtract to u,v to create outliers.  Now multiplied by sd.
#alphaWind: outliers must have errors added in the direction of alpha.  However, the angle
#  need not be alpha exactly, but instead alpha +/- alphaWind/2.  Defaults to 80 degrees.  
m2AngleRestrict=function(n,p,xi,omega,alpha,nu,p.out,k,alphaWind=2*pi/4.5){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  #p.out=p.out/2
  #C=rbinom(2*n,1,p.out)
  C=rep(0,n)
  C[1:(p.out*n)] = 1
  C = sample(C)
  C = c(C, C)
  
  s=1
  cs.m=matrix(C*s,p,2*n,byrow=T)
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  length = NULL; angle = NULL
  for(i in 1:p){
    #Which values correspond to one pressure level?
    filt = (2*i-1):2*i
    length = cbind(length, sample(size=n, sqrt(diag(omega))[filt], replace=T) )
    alphaAngle = atan2(alpha[filt][1], alpha[filt][2])
    angle = cbind(angle, runif(n, alphaAngle-alphaWind/2, alphaAngle+alphaWind/2 ) )
  }
  length = matrix(length, ncol=p)
  angle = matrix(angle, ncol=p)
  k = t(rbind(cos(angle)*length, sin(angle)*length))*k
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#outlying at all higher levels
m3=function(n,p,xi,omega,alpha,nu,p.out,k){
  uv=rmst(n,xi,omega, alpha, nu)
  #u component: p by n
  u=t(uv[,1:p])
  #v component: p by n
  v=t(uv[,(p+1):(2*p)])
  p.out=p.out/2
  
  C=rbinom(2*n,1,p.out)
  
  cout=which(C==1)
  nout=sum(C)
  s=2*rbinom(nout,1,0.5)-1
  cs.m=matrix(0,p,2*n,byrow=T)
  ti=runif(nout,1,p)
  
  for(j in 1:nout){
    part=ti[j]:p
    cs.m[part,cout[j]]=s[j]
  }
  trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])
  
  y=cbind(u,v)+k*cs.m
  list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}









#############################################################################################
#
#
# fastTLE implementations
#
#
#############################################################################################


# data = rnorm(1000)
# data[1:100] = rnorm(100, sd=100)
# qplot(data)
# mean(data); sd(data)
# temp = fast.TLE.normal(data, k=.95)
# temp = fast.TLE.normal(data, k=.9)
# temp = fast.TLE.normal(data, k=.8)
# temp = fast.TLE.normal(data, k=.7)

# data = rst(100, xi=0, omega=1, alpha=6, nu=12)
# qplot(data)
# data[1:10] = rnorm(10, sd=100)
# qplot(data)
# mean(data); sd(data)
# robustST(x=matrix(1,nrow=length(data)), y=data, robust=F)
# robustST(x=matrix(1,nrow=length(data)), y=data, robust=T)
# temp = fast.TLE.ST(data, k=.95)
# temp = fast.TLE.ST(data, k=.9)
# temp = fast.TLE.ST(data, k=.8)
# temp = fast.TLE.ST(data, k=.7)

# data = rmst(100, xi=c(0,0), Omega=diag(c(3,1)), alpha=c(6,-3), nu=12)
# qplot(data[,1], data[,2])
# data[1:3,] = rmst(3, Omega=diag(c(100,100)), alpha=c(6,-3))
# qplot(data[,1], data[,2])
# robustST(x=matrix(1,nrow=NROW(data)), y=data, robust=F)
# robustST(x=matrix(1,nrow=NROW(data)), y=data, robust=T)
# temp = fast.TLE.BST(data, k=.95)
# temp = fast.TLE.BST(data, k=.9)
# temp = fast.TLE.BST(data, k=.8)
# temp = fast.TLE.BST(data, k=.7)