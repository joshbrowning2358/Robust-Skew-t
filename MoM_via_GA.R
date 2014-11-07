source_github("https://raw.githubusercontent.com/rockclimber112358/Random_Code/master/geneticAlgorithm.R")

creater = function(){
  list( xi=runif(1,-10,10)
       ,omega=10^runif(1,-2,2)
       ,alpha=runif(1,-5,5)
       ,nu=10^runif(1,log(4)/log(10),3)
  )
}

#Assumes m1, m2, m3, m4 are all defined
fitness = function(listEl){
  delta = listEl$alpha/sqrt(1+listEl$alpha^2)
  mu = delta*(lgamma(listEl$nu/2-1/2)-lgamma(listEl$nu/2))*sqrt(listEl$nu/pi)
  
  m1Est = listEl$omega*mu + listEl$xi
  m2Est = listEl$omega^2*listEl$nu/(listEl$nu - 2) + 2*listEl$xi*m1Est - listEl$xi^2
  m3Est = listEl$omega^3*mu*(3 - delta^2)*listEl$nu/(listEl$nu - 3) +
    3*listEl$xi*m2Est - 3*listEl$xi^2*m1Est + listEl$xi^3
  m4Est = 3*listEl$omega^4*listEl$nu^2/((listEl$nu - 2)*(listEl$nu - 4)) +
    4*listEl$xi*m3Est - 6*listEl$xi^2*m2Est + 4*listEl$xi^3*m1Est - listEl$xi^4
  out = (1-(m1Est-m1)^2/m1^2) + (m2^2-(m2Est-m2)^2) + (m3^2-(m3Est-m3)^2) + (1-(m4Est-m4)^2/m4^2)
  return( max(0,out) )
}

crossover = function(listEl1, listEl2){
  out = list( xi=sample(c(listEl1$xi, listEl2$xi), size=1)
             ,omega=sample(c(listEl1$omega, listEl2$omega), size=1)
             ,alpha=sample(c(listEl1$alpha, listEl2$alpha), size=1)
             ,nu=sample(c(listEl1$nu, listEl2$nu), size=1) )
  return(out)
}

check = function(listEl){
  listEl$nu>=4 & listEl$omega>0
}

`

mutation = function(listEl){
  mutate = sample(c(0,0,1,1))
  for(i in c(1,3))
    if(mutate[i])
      listEl[i] = listEl[[i]] + rnorm(1,sd=.5)
  for(i in c(2,4))
    if(mutate[i])
      listEl[i] = listEl[[i]]*rnorm(1,mean=1,sd=.2)
  if(listEl[4]<4) listEl[4] = 4
  return(listEl)
}

n=1
xi = runif(n,-10,10)
omega = 10^runif(n,-2,2)
alpha = runif(n,-5,5)
nu = 10^runif(n,log(4)/log(10),3)
delta = alpha/sqrt(1+alpha^2)
mu = delta*(lgamma(nu/2-1/2)-lgamma(nu/2))*sqrt(nu/pi)

m1 = omega*mu + xi
m2 = omega^2*nu/(nu - 2) + 2*xi*m1 - xi^2
m3 = omega^3*mu*(3 - delta^2)*nu/(nu - 3) + 3*xi*m2 - 3*xi^2*m1 + xi^3
m4 = 3*omega^4*nu^2/((nu - 2)*(nu - 4)) + 4*xi*m3 - 6*xi^2*m2 + 4*xi^3*m1 - xi^4

GA( creater, fitness, crossover, mutation=mutation, check=check, popSize=1000, nKeepTop=0)
