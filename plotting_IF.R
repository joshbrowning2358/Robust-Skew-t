library(reshape)
library(gridExtra)
library(numDeriv)
library(sn)
rho = function(x, xi=xi0, omega=omega0, alpha=alpha0, nu=nu0){-log( dst(x, xi, omega, alpha, nu) )}

plotIF = function(xi0=0, omega0=1, alpha0=0, nu0=10000, xSeq = seq(-5,5,.1)){
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
  toPlot = melt(out, id.vars="x")
  toPlot$variable = gsub("psi_", "", toPlot$variable)
  toPlot$variable = factor(toPlot$variable, levels=c("xi", "omega", "alpha", "nu") )
  IF = ggplot( toPlot, aes(x=x, y=value, color=variable) ) + geom_line() +
      labs(y="Psi (proportional to influence function)", title=paste("xi:", xi0, " omega:", omega0, " alpha:", alpha0, " nu:", nu0) )
  return(out)
}

plotIFrobust = function(xi0=0, omega0=1, alpha0=0, nu0=10000, xSeq = seq(-5,5,.1), k=2){
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
  out = data.frame( x=xSeq, psi_xi, psi_omega, psi_alpha, psi_nu )
  toPlot = melt(out, id.vars="x")
  toPlot$variable = gsub("psi_", "", toPlot$variable)
  toPlot$variable = factor(toPlot$variable, levels=c("xi", "omega", "alpha", "nu") )
  IF = ggplot( toPlot, aes(x=x, y=value, color=variable) ) + geom_line() +
      labs(y="Psi (proportional to influence function)", title=paste("xi:", xi0, " omega:", omega0, " alpha:", alpha0, " nu:", nu0) )
  return(out)
}

IF1 = plotIF( xi=0, omega=1, alpha=2, nu=10, xSeq=seq(-20,20,.1) )
p.xi = qplot(IF1$x, IF1$psi_xi, geom="line") + labs(x="y", y="IF(y;xi)", linetype="")
p.omega = qplot(IF1$x, IF1$psi_omega, geom="line") + labs(x="y", y="IF(y;omega)", linetype="")
p.alpha = qplot(IF1$x, IF1$psi_alpha, geom="line") + labs(x="y", y="IF(y;alpha)", linetype="")
p.nu = qplot(IF1$x, IF1$psi_nu, geom="line") + labs(x="y", y="IF(y;nu)", linetype="")
png("~/Professional Files/Mines/Research/Wind Mandy Ying/Documents/st_influence_functions.png", width=8, height=8, units="in", res=600)
grid.arrange( p.xi, p.omega, p.alpha, p.nu, nrow=2, main="Influence functions for ST(0,1,2,10)" )
dev.off()

IF2 = plotIFrobust( xi=0, omega=1, alpha=2, nu=10, xSeq=seq(-20,20,.1), k=8 )
p.xi = qplot(IF1$x, IF1$psi_xi, geom="line", linetype="MLE") + geom_line(data=data.frame(IF2), aes(y=psi_xi, linetype="robust")) + labs(x="y", y="IF(y;xi)", linetype="")
p.omega = qplot(IF1$x, IF1$psi_omega, geom="line", linetype="MLE") + geom_line(data=data.frame(IF2), aes(y=psi_omega, linetype="robust")) + labs(x="y", y="IF(y;omega)", linetype="")
p.alpha = qplot(IF1$x, IF1$psi_alpha, geom="line", linetype="MLE") + geom_line(data=data.frame(IF2), aes(y=psi_alpha, linetype="robust")) + labs(x="y", y="IF(y;alpha)", linetype="")
p.nu = qplot(IF1$x, IF1$psi_nu, geom="line", linetype="MLE") + geom_line(data=data.frame(IF2), aes(y=psi_nu, linetype="robust")) + labs(x="y", y="IF(y;nu)", linetype="")
png("~/Professional Files/Mines/Research/Wind Mandy Ying/Documents/st_robust_influence_functions.png", width=8*12/10, height=8, units="in", res=600)
grid.arrange( p.xi, p.omega, p.alpha, p.nu, nrow=2, main="Influence functions for ST(0,1,2,10)" )
dev.off()
