library(sn)
source("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/Code/sn-funct.R")
source("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/Code/modified_mst.pdev.grad.R")

x = matrix(1,nrow=100)
y = rnorm(100)
dp = c(0,1,0,1000)
llvec = st.pdev(dp, x=x, y=y)
llobs = sapply( y, st.pdev, dp=dp, x=matrix(1) )
llvec
sum( llobs )

dllvec = st.pdev.gh(dp, x=x, y)
dllobs = sapply( y, st.pdev.gh, dp=dp, x=matrix(1) )
dllvec
rowSums(dllobs)

y = matrix(rnorm(200), ncol=2)
param = list( xi=rep(0,NCOL(y)), Omega=diag(NCOL(y)), alpha=rep(0,NCOL(y)), nu=1 )
param = dplist2optpar(param)
llvec = mst.pdev(param, x=x, y=y, w=1)
llobs = c()
for(i in 1:NROW(y))
  llobs = c(llobs, mst.pdev(param, x=matrix(1), y=y[i,,drop=F], w=1))
llvec
sum( llobs )

dllvec = mst.pdev.grad( param, x=x, y=y, w=1)
#Fix: w must be a numeric vector
#dllvec = mst.pdev.grad( param, x=x, y=y, w=rep(1,NROW(y)))
dllobs = matrix(0, nrow=0, ncol=length(dllvec))
for(i in 1:NROW(x))
  dllobs = rbind(dllobs, mst.pdev.grad(param, x=matrix(1), y=y[i,,drop=F], w=1))
dllvec
colSums(dllobs)
dllvec - colSums(dllobs)

#Using new implementation
dllvec = mst.pdev.grad.new( param, x=x, y=y, w=1)
dllobs = matrix(0, nrow=0, ncol=length(dllvec))
for(i in 1:NROW(x))
  dllobs = rbind(dllobs, mst.pdev.grad.new(param, x=matrix(1), y=y[i,,drop=F], w=1))
dllvec
colSums(dllobs)
dllvec - colSums(dllobs)

#How do the mst.pdev.grad functions compare to numerical derivatives?
#param = rep(0,8)
param = rnorm(8)
num.grad = c()
for(i in 1:8)
  num.grad=c(num.grad, grad( function(x){
    param2 = param
    param2[i] = x
    mst.pdev(param2, x=matrix(1,nr=NROW(y)), y=y, w=1)}, x=param[i] ) )

num.grad - mst.pdev.grad(param, x=matrix(1,nr=NROW(y)), y=y, w=1)
num.grad - mst.pdev.grad.new(param, x=matrix(1,nr=NROW(y)), y=y, w=1)

#dp = c(0,1,0,1000)
dp = rnorm(4); dp[c(2,4)] = exp(dp[c(2,4)])
num.grad = c()
for(i in 1:4)
  num.grad=c(num.grad, grad( function(x){
    dp2 = dp
    dp2[i] = x
    st.pdev(dp2, x=matrix(1,nr=NROW(y)), y=y[,1], w=1)}, x=dp[i] ) )

st.pdev.gh( dp, x=matrix(1,nr=NROW(y)), y=y[,1], w=1)
num.grad