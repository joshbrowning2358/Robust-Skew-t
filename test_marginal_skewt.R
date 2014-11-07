setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
source("Code/functions.R")

n = 100000

params = denverParams("EX")
xi = params$xi
omega = params$omega
alpha = params$alpha
nu = params$nu
full = rmst(n, xi=xi, Omega=omega, alpha=alpha, nu=nu)

params2 = marginal(xi, omega, alpha, nu, r=c(1,9), alphaAdj=T )
xi2 = params2$xi
omega2 = params2$omega
alpha2 = params2$alpha
nu2 = params2$nu
margAdj = rmst(n, xi=xi2, Omega=omega2, alpha=alpha2, nu=nu2)

params3 = marginal(xi, omega, alpha, nu, r=c(1,9), alphaAdj=F )
xi3 = params3$xi
omega3 = params3$omega
alpha3 = params3$alpha
nu3 = params3$nu
margNaive = rmst(n, xi=xi3, Omega=omega3, alpha=alpha3, nu=nu3)

# Compare distributions of simulated data:
results = data.frame( full[,c(1,9)], model="16-variate")
results = rbind(results, data.frame(margAdj, model="Marginal"))
results = rbind(results, data.frame(margNaive, model="Naive"))

dir.create("Results/alpha_marginals")
ggsave("Results/alpha_marginals/x1.png",
  ggplot(results, aes(x=X1, color=model) ) + geom_density() + coord_cartesian(xlim=c(-20,20))
)
ggsave("Results/alpha_marginals/x2.png",
  ggplot(results, aes(x=X2, color=model) ) + geom_density() + coord_cartesian(xlim=c(-20,20))
)

ggsave("Results/alpha_marginals/Density_Comparison.png",
  ggplot(results, aes(x=X1, y=X2, color=model) ) + geom_density2d(n=10) +
    facet_wrap( ~ model )
)