##simulation examples
rm(list=ls())
library(fields)
library(mvoutlier)
library(depth)
library(spatialkernel)
library(scatterplot3d)
library(rgl)
source('model.R')
source('method.R')
source('sn/R/sn.R')



##########parameters
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

#Values based on those observed at Denver Station
alpha.obs=c(2.16,  1.44,  1.35,  0.95,  0.60,  0.31,  2.62,  3.02, -0.03, -0.81, -0.82, -0.89, -0.82, -0.75, -0.31,  0.39 )
#alpha.obs=c(0.09, -0.07,  0.13, -0.25,  0.18, -0.08, -0.04,  0.00, -0.09,  0.11,  0.02, -0.08, -0.03,  0.07,  0.07, -0.10)
df.obs=10

#Values for a MVN distribution
alpha.MVN=rep(0,16)
df.MVN=Inf

#Values for distributions that are more skewed than what was observed at the Denver Station
alpha.EX=c(3:10,3:10)
df.EX=5

skew3=cbind(alpha.obs,alpha.MVN,alpha.EX)
df3=c(df.obs,df.MVN,df.EX)
#Sample sizes to loop over
nsim=1000
#n.vals=seq(500,3000,by=500)

p=8

##---------------------------------------------------------

set.seed(10)
#for(i in 1:length(n.vals)){
n=500	
	#n=n.vals[i]
table=NULL
for (mm in 1:3){
tab=NULL
for(j in 1:nsim){
rw=NULL
###########

alpha=skew3[,mm]
df=df3[mm]

#model 2: outlying at all levels
#contamination magnitude 

k=6
std=sqrt(diag(omega2))
kmat=cbind(matrix(std[1:p],p,n),matrix(std[(p+1):(2*p)],p,n))
kk=k*kmat
#outlier percentage 
p.out=0.1

uv=m3(n,p,xi,omega2,alpha,df,p.out,kk)
u=uv$u
v=uv$v
out=uv$outlier

#biv norm
temp=tpfp(mvns(u,v),out)
rw=c(rw,temp)
#mvn
temp=tpfp(mvnf(u,v),out)
rw=c(rw,temp)
#depth
out.dp=dp(u,v,3)
temp=tpfp(out.dp,out)
rw=c(rw,temp)
#biv skew
temp=tpfp(msts(u,v),out)
rw=c(rw,temp)
#mvst
temp=tpfp(mstf(u,v),out)
rw=c(rw,temp)
tab=rbind(tab,rw)
}
table=cbind(table,rbind(apply(tab,2,mean),apply(tab,2,sd)))
}

full=signif(table*100,3)
full[,1:10]
full[,11:20]
full[,21:30]


###################################################
##---------------------------------------------
#Plots of 100 launches of each type
##---------------------------------------------
###################################################

set.seed(100)
mm=2
n=500
alpha=skew3[,mm]
df=df3[mm]


k=6
std=sqrt(diag(omega2))
kmat=cbind(matrix(std[1:p],p,n),matrix(std[(p+1):(2*p)],p,n))
kk=k*kmat
#outlier percentage 
p.out=0.1

uv=m3(n,p,xi,omega2,alpha,df,p.out,kk)
u=uv$u
v=uv$v
out=uv$outlier


start=1
nlaunch=100
colors=rainbow(nlaunch+1)
index=start:(start+nlaunch-1)
z=c(1,3,4,5,5.5,6,7,7.3)

index=start:nlaunch
umax=max(apply(u[,index],2,max),na.rm=TRUE)
umin=min(apply(u[,index],2,min),na.rm=TRUE)

vmax=max(apply(v[,index],2,max),na.rm=TRUE)
vmin=min(apply(v[,index],2,min),na.rm=TRUE)

launch.out=apply(out,2,sum)

plot3d(5,2,1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), type="n",xlab="U",ylab="V",zlab="Pressure Level",main="",top=TRUE,label.tick.marks=FALSE)

for(i in 1:nlaunch){
	if(launch.out[i]==0){
		plot3d(x=u[,i], y=v[,i] , z , col=1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
	if(launch.out[i]>0){
		plot3d(x=u[,i], y=v[,i] , z , col=4, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
		
}
	

decorate3d(main="MVN: Some Contaminated",xlab="",ylab="",cex.main=1.5)


###################################################

set.seed(100)
mm=1
n=500
alpha=skew3[,mm]
df=df3[mm]


k=6
std=sqrt(diag(omega2))
kmat=cbind(matrix(std[1:p],p,n),matrix(std[(p+1):(2*p)],p,n))
kk=k*kmat
#outlier percentage 
p.out=0.1

uv=m3(n,p,xi,omega2,alpha,df,p.out,kk)
u=uv$u
v=uv$v
out=uv$outlier


start=1
nlaunch=100
colors=rainbow(nlaunch+1)
index=start:(start+nlaunch-1)
z=c(1,3,4,5,5.5,6,7,7.3)

index=start:nlaunch
umax=max(apply(u[,index],2,max),na.rm=TRUE)
umin=min(apply(u[,index],2,min),na.rm=TRUE)

vmax=max(apply(v[,index],2,max),na.rm=TRUE)
vmin=min(apply(v[,index],2,min),na.rm=TRUE)

launch.out=apply(out,2,sum)

plot3d(5,2,1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), type="n",xlab="U",ylab="V",zlab="Pressure Level",main="",top=TRUE,label.tick.marks=FALSE)

for(i in 1:nlaunch){
	if(launch.out[i]==0){
		plot3d(x=u[,i], y=v[,i] , z , col=1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
	if(launch.out[i]>0){
		plot3d(x=u[,i], y=v[,i] , z , col=4, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
		
}
	

decorate3d(main="OBS: Some Contaminated",xlab="",ylab="",cex.main=1.5)



###################################################

set.seed(100)
mm=3
n=500
alpha=skew3[,mm]
df=df3[mm]


k=6
std=sqrt(diag(omega2))
kmat=cbind(matrix(std[1:p],p,n),matrix(std[(p+1):(2*p)],p,n))
kk=k*kmat
#outlier percentage 
p.out=0.1

uv=m3(n,p,xi,omega2,alpha,df,p.out,kk)
u=uv$u
v=uv$v
out=uv$outlier


start=1
nlaunch=100
colors=rainbow(nlaunch+1)
index=start:(start+nlaunch-1)
z=c(1,3,4,5,5.5,6,7,7.3)

index=start:nlaunch
umax=max(apply(u[,index],2,max),na.rm=TRUE)
umin=min(apply(u[,index],2,min),na.rm=TRUE)

vmax=max(apply(v[,index],2,max),na.rm=TRUE)
vmin=min(apply(v[,index],2,min),na.rm=TRUE)

launch.out=apply(out,2,sum)

plot3d(5,2,1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), type="n",xlab="U",ylab="V",zlab="Pressure Level",main="",top=TRUE,label.tick.marks=FALSE)

for(i in 1:nlaunch){
	if(launch.out[i]==0){
		plot3d(x=u[,i], y=v[,i] , z , col=1, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
	if(launch.out[i]>0){
		plot3d(x=u[,i], y=v[,i] , z , col=4, xlim=c(umin,umax), ylim=c(vmin,vmax), zlim=c(1,8), xlab="U",ylab="V",zlab="Pressure Level",main="",type="l",add=TRUE,size=10)}
		
}
	

decorate3d(main="EX: Some Contaminated",xlab="",ylab="",cex.main=1.5)





