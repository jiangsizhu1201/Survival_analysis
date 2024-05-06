rm(list=ls())
library(survival)
library(survminer) 
library(ggstatsplot) 

load('cosg_celltype__marker_cosg.Rdata') 
head(marker_cosg)
symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
names(symbols_list)

load( file = 'batch_cox_results.Rdata')
library(remotes)
cox_results = na.omit(cox_results)
do.call(rbind,lapply(symbols_list, function(x){
  # x = symbols_list[[1]]
  x = x[x%in% rownames(cox_results)]
  this_cox = as.data.frame(cox_results[x,])
  table(this_cox$p < 0.01)
}))





