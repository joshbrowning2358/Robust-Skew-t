
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

#outlying at random l levels
m4=function(n,p,xi,omega,alpha,nu,p.out,k,l){
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
for(j in 1:nout){
#random l level
part=sample(1:p,l)
cs.m[part,cout[j]]=s[j]
}
trueout=abs(cs.m[,1:n])|abs(cs.m[,(n+1):(2*n)])

y=cbind(u,v)+k*cs.m
list(u=y[,1:n],v=y[,(n+1):(2*n)],outlier=trueout)
}

#flip
m5=function(n,p,xi,omega,alpha,nu,p.out){
uv=rmst(n,xi,omega, alpha, nu)
#u component: p by n
u=t(uv[,1:p])
#v component: p by n
v=t(uv[,(p+1):(2*p)])
p.out=p.out

C=rbinom(n,1,p.out)
cout=which(C==1)
nout=sum(C)

trueout=matrix(0,p,n)
for(j in 1:nout){

#part=sample(1:p,2)
#temp=u[part[1],cout[j]]
#u[part[1],cout[j]]=u[part[2],cout[j]]
#u[part[2],cout[j]]=temp
#temp=v[part[1],cout[j]]
#v[part[1],cout[j]]=v[part[2],cout[j]]
#v[part[2],cout[j]]=temp
u[,cout[j]]=rev(u[,cout[j]])
v[,cout[j]]=rev(v[,cout[j]])
trueout[,cout[j]]=1
}

list(u=u,v=v,outlier=trueout)
}


#check detected outlying points (u,v)

#for model 1 without outliers
#out: detected outliers; p by n matrix
#n: total number of launches
tp=function(out){
np=prod(dim(out))
if (sum(out)==0) {pright=1;pwrong=0}
else {pright=0;pwrong=sum(out)/(np)}
return(c(pright,pwrong))
}

#for outlier models
#trueout: p by n T/F matrix
tpfp=function(out,trueout){
np=prod(dim(out))
pright=sum(out & trueout)/sum(trueout)
pwrong=sum((out-trueout)==1)/(np-sum(trueout))
return(c(pright,pwrong))
}


#n: total number of launches
tpf=function(out){
outf=colSums(out)>0
if (sum(outf)==0) {pright=1;pwrong=0}
else {pright=0;pwrong=sum(out)/n}
return(c(pright,pwrong))
}

#for outlier models
#trueout: p by n T/F matrix
tpfpf=function(out,trueout){
outf=colSums(out)>0
trueoutf=colSums(trueout)>0
pright=sum(outf & trueoutf)/sum(trueoutf)
pwrong=sum((outf-trueoutf)==1)/(n-sum(trueoutf))
return(c(pright,pwrong))
}















