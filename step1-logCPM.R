rm(list=ls())

library(data.table) 
dat <- data.table::fread('TCGA-LUAD.htseq_counts.tsv.gz',
                          data.table = F) #全部是symbol
head(dat[,1:4])
tail(dat[,1:4]) 
dat = dat[1:(nrow(dat)-5),]
rownames(dat) = dat$Ensembl_ID
a = dat
a = a[,-1]
##逆转 log
a = as.matrix(2^a - 1)
# 用apply转换为整数矩阵
head(a[,1:4])
tail(a[,1:4]) 
colSums(a)/1e6
# 普通转录组定量 20m
# 想看融合基因，可变剪切，100M
exp = apply(a, 2, as.integer)
rownames(exp) = rownames(dat)
exp= log(edgeR::cpm(exp)+1)
library(stringr)
head(rownames(exp))
library(AnnoProbe)
library(tinyarray)
rownames(exp) = substr(rownames(exp), 1, 15) 
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
head(exp[,1:4])
tail(exp[,1:4]) 
proj='tcga-luad'
save(exp,file = paste0(proj,".htseq_counts.rdata") )



 

