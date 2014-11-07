cArgs = commandArgs(trailingOnly=TRUE)
#cArgs = c(1,10)
if(length(cArgs)!=2)
  stop("Exactly two args are required, alpha and nu!")

alpha = as.numeric(cArgs[1])
nu = as.numeric(cArgs[2])
cat("Using alpha=",alpha,"\n")
cat("Using nu=",nu,"\n")

library(sn)
library(ggplot2)
library(scales)
library(plyr)
library(reshape)

if(Sys.info()[1]=="Windows")
  setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Robust Estimators/")
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="jb")
  setwd("/media/storage/Professional Files/Mines/Research/Robust Estimators/")
if(Sys.info()[1]=="Linux" & grepl("ch120", Sys.info()[4]))
  setwd("~/Research/Robust Estimators")
source("Code/sn-funct.R")
source("Code/functions.R")

prefix = paste0("Results/simulation_results_a",alpha,"_n",nu)

##########################################################################
# Run the simulation
##########################################################################

params = data.frame(outPct=c(0,.01,.05))
params = merge(params, data.frame(n=c(30,100,300)))
params = merge(params, data.frame(repl=1:100))
results = matrix(0, nrow=0, ncol=14)
colnames(results) = c("n", "k", "outCnt", "xi", "omega", "alpha", "nu", "xiEst", "omegaEst"
      ,"alphaEst", "nuEst", "compTime", "estimator", "runNo")
for(i in 1:nrow(params)){
  out = runSim(n=params[i,"n"], outPct=params[i,"outPct"], k=4:10, outSigma=5, xi=0, omega=1, alpha=alpha, nu=nu)
  out$runNo = i
  colnames(out)=colnames(results)
  results = rbind(results, out)
  if(i%%10==0 ){
    cat("Run",i,"completed.\n")
    write.csv(results, file=paste0(prefix,".csv"), row.names=F)
  }
}

fname = paste0(prefix,".RData")
save(results, params, file=fname)

##################################################################
# Analyzing results
##################################################################

table( results$estimator[!is.na(results$xiEst)] )/max(results$runNo)
#results = results[results$estimator!="nlmRob",]
noNA = ddply(results, "runNo", function(df){
  !any(is.na(df))
} )
results = merge(results, noNA, by="runNo")
colnames(results)[colnames(results)=="V1"] = "noNA"
noNA2 = ddply(results, "runNo", function(df){
  !any(is.na(df[df$estimator!="nlmRob",]))
} )
results = merge(results, noNA2, by="runNo")
colnames(results)[colnames(results)=="V1"] = "noNA2"

results$outPct = results$outCnt / results$n
results$n = as.numeric(results$n)

ggsave(paste0(prefix,"_convergence.png"),
  ggplot(results, aes(x=estimator, y=runNo, fill=!is.na(xiEst)) ) + geom_tile() +
    labs(fill="Converged?") + scale_y_continuous(limits=c(1,nrow(params)))
  ,width=6, height=8)

plot_results( results[results$noNA2 & results$estimator!="nlmRob",], prefix=paste0(prefix,"_no_NlmRob"), alpha=alpha, nu=nu)
plot_results( results[results$noNA,], prefix=paste0(prefix,"_all"), alpha=alpha, nu=nu)
