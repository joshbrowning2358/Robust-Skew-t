library(stringr)
library(reshape)
library(ggplot2)
setwd("~/Professional Files/Mines/Research/Robust Estimators/")
files = list.files("Results/Simulation 20141023/")
files = paste0("Results/Simulation 20141023/", files[grepl("RData",files)])

table=NULL
tablesd=NULL
for(file in files){
  load(file)
  model = as.numeric(gsub("mod_", "", str_extract( file, "mod_[1-5]" ) ))
  method = as.numeric(gsub("mm_", "", str_extract( file, "mm_[1-3]" ) ))
  outlier = as.numeric(gsub("pout_", "", str_extract( file, "pout_[0-9\\.]{4}" ) ))
  k = as.numeric(gsub("k_", "", str_extract( file, "k_[0-9]{1,2}" ) ))
  if(is.na(k))
    k = 10
  if(is.na(outlier))
    outlier = 0.05
  percent=c(model, method, outlier, k, apply(tab,2,mean), mean(apply(conv, 2, mean)))
  percentsd=c(model, method, outlier, k, apply(tab,2,sd),NA)
  table=rbind(table,percent)
  tablesd=rbind(tablesd,percentsd)
}
colnames(table) = c("mod", "dist", "outliers", "k", "TP mean", "FP mean", "Conv")
table = data.frame(table)
table[,2] = ifelse(table[,2]==1, "MVN", ifelse(table[,2]==2, "OBS", "EX") )
table = cbind(table, tablesd[,5:6])
colnames(table) = c("mod", "dist", "outliers", "k", "TP mean", "FP mean", "Conv", "TP std", "FP std")
table$mod = paste0("Model ", table$mod)
table[,5:9] = table[,5:9]*100
table = melt(table, id.vars=c("mod", "dist", "outliers", "k"))
table = cast(table, "mod + dist + outliers + variable ~ k" )
table$outliers = table$outliers*100
colnames(table)[colnames(table)=="variable"] = "type"

results = read.csv("Results/Simulation 20140903/results_from_ying.txt", header=F
        ,stringsAsFactors=F)
outliers = 10
rm(results2)
for(i in 1:nrow(results)){
#for(i in 1:20){
  if(grepl("Model", results[i,]))
    mod = results[i,]
  if(gsub(" ","",results[i,]) %in% c("MVN", "OBS", "EX"))
    dist = gsub(" ","",results[i,])
  if(grepl("outliers",results[i,]))
    outliers = str_extract(results[i,], "[0-9]*")
  if(grepl("(TP|FP)", results[i,])){
    type = str_extract(results[i,], "^[A-Za-z ]*")
    vals = str_extract(results[i,], "[0-9\\. ]*$")
    vals = strsplit(vals, " ")[[1]]
    vals = as.numeric(vals[vals!=""])
    row = data.frame(mod, dist, outliers, type, vals[1], vals[2], vals[3], vals[4], vals[5])
    if(exists("results2"))
      results2 = rbind(results2, row)
    else
      results2 = row
  }
}
results2$type = gsub(" *$", "", results2$type)

finalResults = merge( results2, table, by=c("mod", "dist", "outliers", "type"), all.y=T )
#finalResults = table #Use if not tying in Ying's output

toPlot = melt( finalResults, id.vars=c("mod", "dist", "outliers", "type") )
toPlot = cast(toPlot, formula=mod + dist + outliers + variable ~ type, value="value"
    ,fun.aggregate=function(x){mean(x)} )
toPlot$variable = as.numeric(gsub("(vals\\.|\\.)", "", toPlot$variable) )
colnames(toPlot) = gsub(" ", "\\.", colnames(toPlot) )
ggplot(toPlot, aes(x=factor(variable), y=TP.mean, color=dist, group=paste0(dist,variable)) ) +
  geom_boxplot() +
  facet_grid( outliers ~ mod )

ggplot(toPlot, aes(x=variable, y=FP.mean, color=dist, group=paste0(dist,variable)) ) +
  geom_boxplot() +
  facet_grid( outliers ~ mod )

ggplot(toPlot, aes(x=variable, y=Conv, color=dist, group=paste0(dist,variable)) ) +
  geom_boxplot() +
  facet_grid( outliers ~ mod )

#Write output in Ying's format:
sink("Results/Simulation 20141023/results_from_josh_k=10_out=5.txt")
for(model in 1:4){
  cat("Model", model, "\n\n")
  for(mm in c("MVN", "OBS", "EX")){
    cat(mm, "\n")
    toPrint = finalResults[finalResults$mod==paste("Model",model) & finalResults$dist==mm
            & finalResults$outliers==5,]
    cat(paste0("TP mean\t", paste(sapply(toPrint[toPrint$type=="TP mean",5:ncol(toPrint)], formatC, digits=3, format="f"),collapse="\t")), "\n")
    cat(paste0("FP mean\t", paste(sapply(toPrint[toPrint$type=="FP mean",5:ncol(toPrint)], formatC, digits=3, format="f"),collapse="\t")), "\n")
    cat("\n")
    cat(paste0("TP std\t", paste(sapply(toPrint[toPrint$type=="TP std",5:ncol(toPrint)], formatC, digits=3, format="f"),collapse="\t")), "\n")
    cat(paste0("FP std\t", paste(sapply(toPrint[toPrint$type=="FP std",5:ncol(toPrint)], formatC, digits=3, format="f"),collapse="\t")), "\n")
    cat("\n")
    cat(paste0("Conv %\t", paste(sapply(100-toPrint[toPrint$type=="Conv",5:ncol(toPrint)], formatC, digits=3, format="f"),collapse="\t")), "\n")
    cat("\n")
  }
  cat("\n")
}
sink()
