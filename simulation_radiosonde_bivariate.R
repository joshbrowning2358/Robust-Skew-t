cArgs = commandArgs(trailingOnly=TRUE)
#cArgs = c("obs", "300", "T")
if(length(cArgs)!=3)
  stop("Exactly two args are required, obs/MVN/EX and pressure level!")
if(!cArgs[1] %in% c("obs", "MVN", "EX"))
  stop("First argument must be obs or MVN or EX!")
if(!cArgs[2] %in% c(700, 500, 400, 300, 250, 200, 100, 70))
  stop("Second argument must be in c(700, 500, 400, 300, 250, 200, 100, 70)!")
if(!cArgs[3] %in% c("T","F"))
  stop("Third argument must be in c(T,F)!")

type = cArgs[1]
pressure = as.numeric(cArgs[2])
restrict = as.logical(cArgs[3])
cat("Using type=",type,"\n")
cat("Using pressure=",pressure,"\n")
cat("Using restricted angle=",restrict,"\n")

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
source("Code/functions.R")

prefix = paste0("Results/Simulation 20141107/type_", type, "_pressure_",pressure,"_outlier_restrict_",restrict)

##########################################################################
# Run the simulation
##########################################################################

params = data.frame(outPct=c(0,.01,.05))
params = merge(params, data.frame(n=c(100,300,500)))
params = merge(params, data.frame(repl=1:200))
for(i in 1:nrow(params)){
  out = runSim(n=params[i,"n"], outPct=params[i,"outPct"], k=seq(4,26,2), outSigma=5
      ,pressure=pressure, type=type, restrict=restrict)
  out$runNo = i
  if(exists("results"))
    results = rbind(results, out)
  else
    results = out

  if(i%%1==0 )
    cat("Run",i,"completed out of",nrow(params),"\n")
  if(i%%10==0 )
    save(results, params, file=paste0(prefix,".RData"))
}

save(results, params, file=paste0(prefix,".RData"))

##################################################################
# Analyzing results
##################################################################
# 
# table( results$estimator[!is.na(results$xiEst)] )/max(results$runNo)
# #results = results[results$estimator!="nlmRob",]
# noNA = ddply(results, "runNo", function(df){
#   !any(is.na(df))
# } )
# results = merge(results, noNA, by="runNo")
# colnames(results)[colnames(results)=="V1"] = "noNA"
# noNA2 = ddply(results, "runNo", function(df){
#   !any(is.na(df[df$estimator!="nlmRob",]))
# } )
# results = merge(results, noNA2, by="runNo")
# colnames(results)[colnames(results)=="V1"] = "noNA2"
# 
# results$outPct = results$outCnt / results$n
# results$n = as.numeric(results$n)
# 
# ggsave(paste0(prefix,"_convergence.png"),
#   ggplot(results, aes(x=estimator, y=runNo, fill=!is.na(xiEst)) ) + geom_tile() +
#     labs(fill="Converged?") + scale_y_continuous(limits=c(1,nrow(params)))
#   ,width=6, height=8)
# 
# plot_results( results[results$noNA2 & results$estimator!="nlmRob",], prefix=paste0(prefix,"_no_NlmRob"), alpha=alpha, nu=nu)
# plot_results( results[results$noNA,], prefix=paste0(prefix,"_all"), alpha=alpha, nu=nu)
