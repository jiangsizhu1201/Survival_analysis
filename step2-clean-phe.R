rm(list=ls())

proj='tcga-luad'
load(file = paste0(proj,".htseq_counts.rdata") )
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')  
Group = factor(Group,levels = c("normal","tumor"))
print(table(Group))
# 生存分析只需要tumor样品即可
exprSet = exp[,Group=='tumor']


clinical = read.delim('TCGA-LUAD.GDC_phenotype.tsv.gz',
                      fill = T,header = T,sep = "\t")
surv = read.delim('TCGA-LUAD.survival.tsv',header = T) 
library(tidyverse)
meta = left_join(surv,clinical,by = c("sample"= "submitter_id.samples"))
head(meta[,1:4])
tail(meta[,1:4]) 
print(dim(meta))

#### 3.3 样本过滤
#去掉生存信息不全或者生存时间小于30天的样本，样本纳排标准不唯一，且差别很大.
k1 = meta$OS.time >= 30
k2 = !(is.na(meta$OS.time)|is.na(meta$OS))
meta = meta[k1&k2,]
meta = meta[,c(
  'sample',
  'OS',
  'OS.time'
)]
colnames(meta)=c('ID','event','time')
meta$time = meta$time/30
rownames(meta) <- meta$ID
s = intersect(rownames(meta),colnames(exprSet))
exprSet = exprSet[,s]
meta = meta[s,]
identical(rownames(meta),colnames(exprSet))
save(exprSet,meta,file = paste0(proj,".for_survival.rdata") )


## meta的行名和exprSet的列名都是样本名，可以对应起来
# 取第一个基因测试看看
library(survival) 
library(survminer)
meta$gene = ifelse(as.numeric( exprSet[1,])> median(as.numeric( exprSet[1,])),'high','low')
# TTTY4C
g = as.numeric( exprSet['TTTY4C',]);g
meta$gene = ifelse( g >   median( g),'high','low')

meta = as.data.frame(meta)
sfit=survfit(Surv(time, event)~gene, data=meta)
ggsurvplot(sfit,
           pval = TRUE)

