setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators")
source("Code/functions.R")

params = denverParams(type="obs")
params700 = marginal( xi=params$xi, omega=params$omega
        , alpha=params$alpha, nu=params$nu, r=c(1,9) )
params300 = marginal( xi=params$xi, omega=params$omega
        , alpha=params$alpha, nu=params$nu, r=c(1,9)+3 )
params100 = marginal( xi=params$xi, omega=params$omega
        , alpha=params$alpha, nu=params$nu, r=c(1,9)+6 )

d = m2(1000, p=1, xi=params700$xi, omega=params700$omega, alpha=params700$alpha
   ,nu=params700$nu, p.out=.05, k=20 )
d = do.call("data.frame", d)
ggsave( "Results/700_outliers.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2(1000, p=1, xi=params300$xi, omega=params300$omega, alpha=params300$alpha
   ,nu=params300$nu, p.out=.05, k=20 )
d = do.call("data.frame", d)
ggsave( "Results/300_outliers.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2(1000, p=1, xi=params100$xi, omega=params100$omega, alpha=params100$alpha
   ,nu=params100$nu, p.out=.05, k=20 )
d = do.call("data.frame", d)
ggsave( "Results/100_outliers.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Adj(1000, p=1, xi=params700$xi, omega=params700$omega, alpha=params700$alpha
   ,nu=params700$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/700_outliers_new.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Adj(1000, p=1, xi=params300$xi, omega=params300$omega, alpha=params300$alpha
   ,nu=params300$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/300_outliers_new.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Adj(1000, p=1, xi=params100$xi, omega=params100$omega, alpha=params100$alpha
   ,nu=params100$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/100_outliers_new.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Angle(1000, p=1, xi=params700$xi, omega=params700$omega, alpha=params700$alpha
   ,nu=params700$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/700_outliers_angle.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Angle(1000, p=1, xi=params300$xi, omega=params300$omega, alpha=params300$alpha
   ,nu=params300$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/300_outliers_angle.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2Angle(1000, p=1, xi=params100$xi, omega=params100$omega, alpha=params100$alpha
   ,nu=params100$nu, p.out=.05, k=10 )
d = do.call("data.frame", d)
ggsave( "Results/100_outliers_angle.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2AngleRestrict(1000, p=1, xi=params700$xi, omega=params700$omega, alpha=params700$alpha
   ,nu=params700$nu, p.out=.05, k=runif(1000,12,16) )
d = do.call("data.frame", d)
ggsave( "Results/700_outliers_angle_restricted.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2AngleRestrict(1000, p=1, xi=params300$xi, omega=params300$omega, alpha=params300$alpha
   ,nu=params300$nu, p.out=.05, k=runif(1000,12,16) )
d = do.call("data.frame", d)
ggsave( "Results/300_outliers_angle_restricted.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )

d = m2AngleRestrict(1000, p=1, xi=params100$xi, omega=params100$omega, alpha=params100$alpha
   ,nu=params100$nu, p.out=.05, k=runif(1000,12,16) )
d = do.call("data.frame", d)
ggsave( "Results/100_outliers_angle_restricted.png", ggplot(d, aes(x=u, y=v, color=outlier) ) + geom_point() )
