#coordinates
cart2polar= function(xy) {
  d =complex(real = xy[,1], imaginary = xy[,2])
  data.frame(size = Mod(d), angle = Arg(d) %% (2*pi))
}

polar2cart=function(xy){
data.frame(x=xy[,1]*cos(xy[,2]),y=xy[,1]*sin(xy[,2]))
}


#sbiv function: compute depth values for each pressure level
#n: number of launches
#uc: u component
#vc: v component
#l: lth pressure level
sbiv=function(u,v,l){
n=ncol(u)
u0=u[l,]
v0=v[l,]
#remove NA
sel=(!is.na(u0))&(!is.na(v0))
if (sum(sel)>0){
uu=u0[sel]
vv=v0[sel]
one=NULL
onefull=rep(NA,n)
for (j in 1:sum(sel)){

one=c(one,depth(c(uu[j],vv[j]),method='Tukey',cbind(uu,vv)))
}
onefull[sel]=one
}
onefull
}


#find outliers given depth values
#dval is a p by n vector. For weighted depth value, dval has p equal rows
#F: factor for outlier detection
#tau: 50% central region
outdp=function(dval,u,v,F,tau=0.5){
p=nrow(u)
n=ncol(u)
outpt=matrix(0,p,n)
for (i in 1:p) {
#order depth values
index=order(dval[i,],decreasing=T)
#find central region for each pressure level 
u0=u[i,]
v0=v[i,]
sel=(!is.na(u0))&(!is.na(v0))
m=ceiling(tau*sum(sel))
partid=intersect(index,which(sel))
center=partid[1:m]
uc=u0[center]
vc=v0[center]
uf=u0[sel]
vf=v0[sel]
cpoint=c(uc[1],vc[1])

bder=chull(uc, vc)
evlp=cbind(uc[bder],vc[bder])
tocenter=t(evlp)-cpoint
pl=cart2polar(t(tocenter))

pl.f=pl[,1]*(1+F)
fen=cpoint+t(polar2cart(cbind(pl.f,pl[,2])))
fence=t(fen)
check=pinpoly(fence,cbind(uf,vf))
idsel=which(sel)
outpt[i,idsel[check==0]]=1
}
return(out=outpt)
}


wbiv=function(dp) {
#number of points on each pressure level 
nplevel=rowSums(!is.na(dp))
wt=nplevel/sum(nplevel)
#weighted average for functional depth
fundp=apply(wt*dp,2,mean,na.rm=T)
fundp
}





#Data depth methods: single level and functional approach
#u,v: p by n matrix, p pressure level; n launches
dp=function(u,v,f1,keep) {
np=dim(u)
p=np[1]
n=np[2]

#compute depth value for each pressure level
sdval=NULL
trim=NULL
for (i in 1:p){
value=sbiv(u,v,i)
sdval=rbind(sdval,value)
m=ceiling(n*keep)
temp=order(value,decreasing=T)
trim=rbind(trim,temp[1:m])
}

#compute depth value as weighted average
fdval=matrix(wbiv(sdval),p,n,byrow=T)
#singel level method
out.s=outdp(sdval,u,v,f1)
#functional approach
#out.f=outdp(fdval,u,v,f2)
return(list(out.s=out.s,trim=trim))
}

#Bivariate normal; single level alpha=0.025
mvns=function(u,v,alpha=0.025){
p=nrow(u)
n=ncol(u)

out=matrix(0,p,n)
for (i in 1:p) {
uv=cbind(u[i,],v[i,])
initial.mean=apply(uv,2,mean)
initial.scatter=cov(uv)
update.stats3=arw(uv,initial.mean,initial.scatter,alpha)
#update.stats3$w is a vector of T/F values.  Flagged observations are assigned a value of "FALSE"
out[i,!update.stats3$w]=1
}
#return outlying launches
return(out=out)
}

#MVN; functional
mvnf=function(u,v,alpha=0.025){
p=nrow(u)
n=ncol(u)
out=matrix(0,p,n)
uv=t(rbind(u,v))

initial.mean=apply(uv,2,mean)
initial.scatter=cov(uv)

update.stats1=arw(uv,initial.mean,initial.scatter,alpha)
m=update.stats1$m
sig=update.stats1$c #16 by 16

for (j in 1:p){
sel=c(j,j+p)
uv=cbind(u[j,],v[j,])
m0=m[sel]
c0=sig[sel,sel]

update.stats=arw(uv,m0,c0,alpha)
#update.stats3$w is a vector of T/F values.  Flagged observations are assigned a value of "FALSE"
out[j,!update.stats$w]=1
}
return(out=out)
}

#MVST M distance

mstpval=function(uv,pars){
n=nrow(uv)
d=ncol(uv) #d=8 or 16
y=uv

X <- rep(1,n) 
X   <- as.matrix(X)
qrX <- qr(X)
beta  <- pars$beta
Omega <- pars$Omega
alpha <- pars$alpha
df    <- pars$df
#df=dftrue
xi    <- X %*% beta
pp  <- d * qf((1:n)/(n+1),d,df)
pp2 <- qchisq((1:n)/(n+1),d)
Xb  <- qr.fitted(qrX,y)
res <- qr.resid(qrX,y)
rad.n  <- apply(res    * (res %*% pd.solve(var(res))), 1, sum)
rad.st <- apply((y-xi) * ((y-xi) %*% pd.solve(Omega)), 1, sum)
prob <- pf(rad.st/d,d,df)
return(list(pval=1-prob,mdist=rad.st))
}

#BST non robust or sweeped

mst=function(u,v,af=0.025,rm.out=F,part){
np=dim(u)
p=np[1]
n=np[2]
uc=NULL
vc=NULL
if (rm.out){
#remove outlying launches
#keep 95%
for (i in 1:p){
uc=rbind(uc,u[i,part[i,]])
vc=rbind(vc,v[i,part[i,]])
}

}else{
uc=u
vc=v

}
out=matrix(0,p,n)

for (j in 1:p){

uvc=cbind(uc[j,],vc[j,])
uv=cbind(u[j,],v[j,])
mle <- mst.mle(y=uvc)
pars=list(beta=mle$dp$beta,Omega=mle$dp$Omega,alpha=mle$dp$alpha,df=mle$dp$df)
pval=mstpval(uv,pars)$pval
out[j,pval<af]=1
}
return(out=out)
}

#BST robust
#Adjusted this function to also return a vector indicating convergence
#...: used to pass additional arguments to robustST()
mst.robust=function(u,v,af=0.025, ...){
  np=dim(u)
  p=np[1]
  n=np[2]
  
  out=matrix(0,p,n)
  
  conv = rep(NA,p)
  for (j in 1:p){
  
    uv=cbind(u[j,],v[j,])
    rest=robustST(uv, ...)
    conv[j] = rest$convergence
    pars=list(beta=rest$beta,alpha=rest$alpha,df=rest$df,Omega=rest$Omega)
    
    pval=mstpval(uv,pars)$pval
    out[j,pval<af]=1
  }
  return(list(out=out, convergence=conv))
}


#BST true

mst.true=function(u,v,af=0.025,pars){
np=dim(u)
p=np[1]
n=np[2]

out=matrix(0,p,n)

for (j in 1:p){

uv=cbind(u[j,],v[j,])
pval=mstpval(uv,pars[[j]])$pval
out[j,pval<af]=1
}
return(out=out)
}
