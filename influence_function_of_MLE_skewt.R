library(sn)
library(ggplot2)
library(reshape)
library(gridExtra)

setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
plotDir = "Results/"
dir.create(plotDir)

xi0 = 2
omega0 = 3
alpha0 = 2
nu0 = 1000

rho = function(x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0){-log( dst(x, xi, omega, alpha, nu) )}
xGrid = seq(-10,20,.1)
ll = qplot( xGrid, rho(xGrid), geom="line" ) + labs(x="x", y="Negative Log-Likelihood" )
dens = qplot( xGrid, dst(xGrid, xi0, omega0, alpha0, nu0), geom="line" ) + labs(x="x", y="Density")
p1 = arrangeGrob(ll, dens, nrow=1, main=paste(" xi:", xi0, " omega:", omega0, " alpha:", alpha0, " nu:", nu0))
ggsave(paste0(plotDir,"skew_t_ll_density.png"), p1)

computeIF = function(xi0=0, omega0=1, alpha0=0, nu0=10000, xSeq = seq(-5,5,.1)){
  psi_xi = sapply( xSeq, function(x){
    -grad( func=function(xi0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=xi0 )
  } )
  psi_omega = sapply( xSeq, function(x){
    -grad( func=function(omega0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=omega0 )
  } )
  psi_alpha = sapply( xSeq, function(x){
    -grad( func=function(alpha0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=alpha0 )
  } )
  psi_nu = sapply( xSeq, function(x){
    -grad( func=function(nu0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=nu0 )
  } )
  out = data.frame( x=xSeq, psi_xi, psi_omega, psi_alpha, psi_nu )
  return(out)
#  out$dens = sapply( xSeq, function(x){dst(x, xi0, omega0, alpha0, nu0)} )
#  dens = ggplot(out, aes(x=x, y=dens) ) + geom_bar(stat="identity") + labs(y="Density")
#  grid.arrange(IF, dens, nrow=1, main=paste("xi:", xi0, " omega:", omega0, " alpha:", alpha0, " nu:", nu0))
}

#rho(x) -> ifelse(rho(x)<k,rho(x),2k-ke^(-rho(x)/k+1))
#psi(x) -> ifelse(rho(x)<k,rho'(x), e^(-rho(x)/k+1) rho'(x)
computeIFrobust = function(xi0=0, omega0=1, alpha0=0, nu0=10000, xSeq = seq(-5,5,.1), k=2){
  psi_xi = sapply( xSeq, function(x){
    ifelse(rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)<k, 1, exp(-rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)/k+1) )*
      -grad( func=function(xi0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=xi0 )
  } )
  psi_omega = sapply( xSeq, function(x){
    ifelse(rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)<k, 1, exp(-rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)/k+1) )*
      -grad( func=function(omega0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=omega0 )
  } )
  psi_alpha = sapply( xSeq, function(x){
    ifelse(rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)<k, 1, exp(-rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)/k+1) )*
      -grad( func=function(alpha0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=alpha0 )
  } )
  psi_nu = sapply( xSeq, function(x){
    ifelse(rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)<k, 1, exp(-rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)/k+1) )*
      -grad( func=function(nu0){rho(x=x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0)}, x=nu0 )
  } )
  out = data.frame( x=xSeq, psi_xi_robust=psi_xi, psi_omega_robust=psi_omega
                  ,psi_alpha_robust=psi_alpha, psi_nu_robust=psi_nu )
  return(out)
}

plotIF = function(xi0=0, omega0=1, alpha0=0, nu0=10000, k=4, ...){
  MLE = computeIF( xi0, omega0, alpha0, nu0, ... )
  MLE$type = "MLE"
  robust = computeIFrobust( xi0, omega0, alpha0, nu0, k, ... )
  robust$type = "Robust"
  colnames(robust) = colnames(MLE)
  toPlot = rbind(MLE, robust)
  toPlot = melt(toPlot, id.vars=c("x","type"))
  toPlot$variable = gsub("psi_", "", toPlot$variable)
  toPlot$variable = factor(toPlot$variable, levels=c("xi","omega","alpha","nu"))
  p = ggplot(toPlot, aes(x=x, linetype=type, color=variable)) + geom_line(aes(y=value)) +
    facet_wrap( ~ variable, scale="free" ) +
    labs(y="Influence function", linetype="Estimator", color="Parameter"
         ,title=paste("xi:", xi0, "omega:", omega0, "alpha:", alpha0, "nu:", nu0))
  return(p)
}

for(k in 0:10)
  ggsave(paste0(plotDir,"Influence_functions_k=",k,".png"),
         plotIF( xi=0, omega=1, alpha=3, nu=3, xSeq=seq(-10,10,.2), k=k )
    ,width=7, height=7)

plotIF( xi=0, omega=1, alpha=0, nu=5, xSeq=seq(-20,20,.2), k=10000 )
plotIF( xi=0, omega=1, alpha=0, nu=5, xSeq=seq(-100,100,.2), k=10000 )
plotIF( xi=0, omega=1, alpha=2, nu=5, xSeq=seq(-20,20,.2), k=10000 )
plotIF( xi=0, omega=1, alpha=2, nu=20, xSeq=seq(-20,20,.2), k=10000 )
plotIF( xi=0, omega=1, alpha=2, nu=100, xSeq=seq(-20,20,.2), k=10000 )
plotIF( xi=0, omega=1, alpha=2, nu=100, xSeq=seq(-100,100,.5), k=10000 )
