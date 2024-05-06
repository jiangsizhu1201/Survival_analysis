rm(list=ls())

library(data.table) 
dat <- data.table::fread('TCGA-LUAD.htseq_counts.tsv.gz',
                          data.table = F) # all symbol

dat = dat[1:(nrow(dat)-5),]
rownames(dat) = dat$Ensembl_ID
a = dat
a = a[,-1]
##transpose log

a = as.matrix(2^a - 1)
# integer matrix

colSums(a)/1e6
# regular transcription quantification 20m
# fusion gene and splice variants 100M

exp = apply(a, 2, as.integer)
rownames(exp) = rownames(dat)
exp= log(edgeR::cpm(exp)+1)

library(stringr)
library(AnnoProbe)
library(tinyarray)

rownames(exp) = substr(rownames(exp), 1, 15) 
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")

proj='tcga-luad'
save(exp,file = paste0(proj,".htseq_counts.rdata") )



 

