library(sn)
library(ggplot2)

?selm.fit #After digging into code, selm.fit calls st.mple/mst.mple to fit parameters
edit(selm.fit) #See line 80 or 157 for call to st.mple/mst.mple
?st.mple #Skew-t maximum penalized likelihood estimation
edit(st.mple)
#st.mple optimizes the log-likelihood via a call to nlminb (the standard optimization function for R).
#To update this with robust estimators, we'll need to modify st.pdev and
#st.pdev.gh, the skew-t penalized deviance and it's gradient.  To robust-ify these
#functions, we'll need to adjust the log-likelihood using a redescending psi function.

#First, try this approach on normal observations
y = rnorm(100)
nlminb(start=c(0,1), function(dp){n.dev(dp, y)}
       ,grad=function(dp){n.dev.gh(dp, y)}
       ,lower=c(-Inf,0))$par
nlminb(start=c(0,1), function(dp){n.dev.robust(dp,y)}
       ,grad=function(dp){n.dev.gh.robust(dp, y)}
       ,lower=c(-Inf,0))$par

y = c(y, 1000)
nlminb(start=c(0,1), function(dp){n.dev(dp, y)}
       ,grad=function(dp){n.dev.gh(dp, y)}
       ,lower=c(-Inf,0))$par
nlminb(start=c(0,1), function(dp){n.dev.robust(dp,y)}
       ,grad=function(dp){n.dev.gh.robust(dp, y)}
       ,lower=c(-Inf,0))$par

####################################################################
# Let's use normal data, but try fitting a skew-t
####################################################################

y = rnorm(1000)
dp = runif(4); dp
nlminb(start=dp, function(dp){st.pdev(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,grad=function(dp){st.pdev.gh(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,lower=c(-Inf,0,-Inf,0))$par
nlminb(start=dp, function(dp){st.pdev.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,grad=function(dp){st.pdev.gh.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,lower=c(-Inf,0,-Inf,0))$par

y = c(y, 1000)
nlminb(start=dp, function(dp){st.pdev(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,grad=function(dp){st.pdev.gh(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,lower=c(-Inf,0,-Inf,0))$par
continue = TRUE
while(continue){
  fit = try(nlminb(start=dp, function(dp){st.pdev.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
         ,grad=function(dp){st.pdev.gh.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
         ,lower=c(-Inf,0,-Inf,0))$par)
  if(!is(fit,"try-error")){
    continue = FALSE
    print(fit)
  }
  dp = runif(4)
}

fit = nlminb(start=dp, function(dp){st.pdev.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,grad=function(dp){st.pdev.gh.robust(dp, y=y, x=matrix(1,nrow=NROW(y)))}
       ,lower=c(-Inf,0,-Inf,0))$par
