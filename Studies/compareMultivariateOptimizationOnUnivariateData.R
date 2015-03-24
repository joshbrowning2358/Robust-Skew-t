## I rewrote the optimization function for fitting the skew-t to just use the
## mst.pdev and mst.pdev.grad function, rather than using one function in the
## univariate case and a different function in the multivariate case.  So, I
## want to compare my optimization with the default sn optimization.
## 
## Note: the default optimization returns parameters of the form
## (xi, omega, alpha, nu)
## and the multivariate optimization returns parameters of the form
## (xi, Omega, alpha, log(nu))
## where Omega is the variance-covariance matrix and omega is the s.d.

myFits = NULL
snFits = NULL
for(i in 1:100){
    y = rnorm(100)
    ## Set pValue really small so that no robustness is included
    myFit = robustSTOnceK(y = y, family = "ST", pValue = .00001,
                          method = "nlminb")
    myFit = c(myFit$xi, myFit$Omega, myFit$alpha, myFit$nu)
    myFit[4] = exp(myFit[4])
    names(myFit) = c("xi", "Omega", "alpha", "nu")
    snFit = sn::selm.fit(x = matrix(1, nrow = NROW(y)), y = y,
                         selm.control = list(), family = "ST")$param$dp
    snFit[2] = snFit[2]^2
    names(snFit)[2] = "Omega"
    myFits = rbind(myFits, myFit)
    snFits = rbind(snFits, snFit)
}

myFits = data.frame(myFits)
snFits = data.frame(snFits)
myFits$type = "new"
snFits$type = "sn"
toPlot = rbind(myFits, snFits)

## Plots seem better in all cases
ggplot(toPlot, aes(x = type, y = xi)) + geom_boxplot()
ggplot(toPlot, aes(x = type, y = Omega)) + geom_boxplot()
ggplot(toPlot, aes(x = type, y = alpha)) + geom_boxplot()
ggplot(toPlot, aes(x = type, y = nu)) + geom_boxplot() +
    scale_y_log10()




## Now let's look at a more complex distribution:
myFits = NULL
snFits = NULL
for(i in 1:100){
    y = sn::rst(n = 300, alpha = 1.3, nu = 25)
    ## Set pValue really small so that no robustness is included
    myFit = robustSTOnceK(y = y, family = "ST", pValue = .00001,
                          method = "nlminb")
    myFit = c(myFit$xi, myFit$Omega, myFit$alpha, myFit$nu)
    myFit[4] = exp(myFit[4])
    names(myFit) = c("xi", "Omega", "alpha", "nu")
    snFit = sn::selm.fit(x = matrix(1, nrow = NROW(y)), y = y,
                         selm.control = list(), family = "ST")$param$dp
    snFit[2] = snFit[2]^2
    names(snFit)[2] = "Omega"
    myFits = rbind(myFits, myFit)
    snFits = rbind(snFits, snFit)
}

myFits = data.frame(myFits)
snFits = data.frame(snFits)
myFits$type = "new"
snFits$type = "sn"
toPlot = rbind(myFits, snFits)

## Plots seem better in all cases
ggplot(toPlot, aes(x = type, y = xi)) + geom_boxplot()
ggplot(toPlot, aes(x = type, y = Omega)) + geom_boxplot()
ggplot(toPlot, aes(x = type, y = alpha)) + geom_boxplot() +
    geom_hline(yintercept = 1.3, color = "red", linetype = 2)
ggplot(toPlot, aes(x = type, y = nu)) + geom_boxplot() +
    scale_y_log10() +
    geom_hline(yintercept = 25, color = "red", linetype = 2)
