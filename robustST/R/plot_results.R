##' Plot results
##' 
##' This function takes the results object as created by runSim and generates
##' plots summarizing the simulation.  These plots are saved in the users
##' current working directory.
##' 
##' @param results A results object, as returned by runSim.
##' @param prefix A prefix to be added to the saved plots.
##' @param xi A parameter of the skew-t distribution.  This should be the
##' parameter used in the simulation.
##' @param omega A parameter of the skew-t distribution.  This should be the
##' parameter used in the simulation.
##' @param alpha A parameter of the skew-t distribution.  This should be the
##' parameter used in the simulation.
##' @param nu A parameter of the skew-t distribution.  This should be the
##' parameter used in the simulation.
##' 
##' @return No results are returned, but plots are saved in the user's current
##' directory.
##' 
##' @export
##' 

plot_results = function(results, prefix, xi=0, omega=1, alpha=0, nu=10000){
    toPlot = results
    toPlot$n = paste0("n=", toPlot$n)
    #Ensure ordering is appropriate
    toPlot$n = factor(toPlot$n, levels=paste0("n=", sort(unique(results$n))) )
    toPlot$outPct = paste0("Outliers: ",toPlot$outPct,"%")
    #Ensure ordering is appropriate
    toPlot$outPct = factor(toPlot$outPct, levels=paste0("Outliers: ", sort(unique(results$outPct)), "%") )
    
    myTheme = theme( axis.text.x=element_text(size=14)
                     ,axis.title.x=element_text(size=14)
                     ,axis.text.y=element_text(size=14)
                     ,axis.title.y=element_text(size=14)
                     ,strip.text.x=element_text(size=14)
                     ,strip.text.y=element_text(size=14))
    
    p = ggplot(toPlot, aes(y=xiEst, x=estimator) ) + geom_boxplot() +
        facet_grid( n ~ outPct ) + geom_hline(yintercept=xi, color="red", linetype=2) +
        myTheme
    ggsave( paste0(prefix, "_xi_est.png"), p, width=8, height=10 )
    
    p = ggplot(toPlot, aes(y=omegaEst, x=estimator) ) + geom_boxplot() +
        facet_grid( n ~ outPct ) + geom_hline(yintercept=omega, color="red", linetype=2) +
        myTheme
    ggsave( paste0(prefix, "_omega_est.png"), p, width=8, height=10 )
    
    p = ggplot(toPlot, aes(y=alphaEst, x=estimator) ) + geom_boxplot() +
        facet_grid( n ~ outPct ) + geom_hline(yintercept=alpha, color="red", linetype=2) +
        myTheme
    ggsave( paste0(prefix, "_alpha_est.png"), p, width=8, height=10 )
    
    p = ggplot(toPlot, aes(y=nuEst, x=estimator) ) + geom_boxplot() +
        facet_grid( n ~ outPct ) + geom_hline(yintercept=nu, color="red", linetype=2) +
        scale_y_log10(breaks=10^(1:10)) + myTheme
    ggsave( paste0(prefix, "_nu_est.png"), p, width=8, height=10 )
}