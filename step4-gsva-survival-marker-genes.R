rm(list=ls())
library(survival)
library(survminer) 
library(ggstatsplot) 
library(gplots)
library(ggplot2) 
library(pheatmap)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)
 
# 1. 载入表达量矩阵和临床信息 ----
proj='tcga-luad'
load(file = paste0(proj,".for_survival.rdata") ) 
exprSet[1:4,1:4]
phe=meta
head(phe)
mySurv <- with(phe, Surv(time, event))
survival_dat=phe

# 2. creat geneset----
# load('../../paper-figures/cosg_celltype__marker_cosg.Rdata')
load('cosg_celltype__marker_cosg.Rdata')
head(marker_cosg$names)
deg_list =  as.list(marker_cosg$names)
names(deg_list)
deg_list
gs = lapply(deg_list, toupper) 
geneset <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
geneset

# 3. run gsva----
X=as.matrix(exprSet)
es.max <- gsva(X, geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz=4)
es.max[1:4, 1:4] 
pheatmap(es.max) 
## gsva clusters
p = pheatmap(es.max,
             show_colnames = F)
p
dev.off()
pbmc_g = cutree(p$tree_col, k=3)
table(pbmc_g)


# 4. 根据gsva结果高低分组后批量生存分析 ----
es.max[1:4, 1:4]
splots <- list()
g = 1
for (i in  names(deg_list) ) {
  # i =  names(deg_list) [1]
  subset = paste0('cluster_',i)
  print(subset)
  v = as.numeric(es.max[i,])   #每一个亚群表达量。
  sub_group <- ifelse( v < 0,"low","high")   #如果表达量小于0的话，就定义为low。gsva处理过表达量。0.几左右
  table(sub_group) 
  phe$sub_group=sub_group
  # Fit survival curves
  require("survival")
  fit <- survfit(Surv(time, event) ~ sub_group, data = phe)
  library("survminer")
  survp <- ggsurvplot(fit, data = phe,
                      surv.median.line = "hv", # Add medians survival
                      pval = TRUE,             # Add p-value and tervals 
                      conf.int = TRUE,        # Add the 95% confidence band
                      risk.table = TRUE,      # Add risk table
                      tables.height = 0.2,
                      tables.theme = theme_cleantable(),
                      palette = "jco",
                      ggtheme = theme_bw(),
                      title = subset)
  print(survp)
  splots[[g]] <-  survp
  g = g + 1
}

length(splots)
x1 = ceiling(sqrt(length(splots)))
y1 = x1

all_plot <- arrange_ggsurvplots(splots,
                                print = F,
                                ncol = x1, nrow = y1,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
# all_plot 
x2=5*x1
y2=5*y1
prefix=''
pro=''
ggsave(all_plot, #path = prefix,
       filename = paste0(pro, 'all_survival_plot.pdf'),
       width = x2,height = y2)

## cut point
## cutpoint
head(phe)
csplots <- list()
cg = 1
for (i in  names(deg_list) ) {
  # i =  names(deg_list) [1]
  subset = paste0('cluster_',i)
  print(subset)
  v = as.numeric(es.max[i,])   #每一个亚群表达量。
  phe$v <- v
  head(phe)
  sur.cut <- surv_cutpoint(phe,
                           time= 'time',
                           event = 'event' ,
                           variables = 'v' )
  sur.cat <- surv_categorize(sur.cut)
  head(sur.cat)
  sfit <- survfit(Surv(time, event) ~ v, data = sur.cat)
  p_surv_cut <- ggsurvplot(sfit, data = phe,
                           surv.median.line = "hv", # Add medians survival
                           pval = TRUE,             # Add p-value and tervals 
                           conf.int = TRUE,        # Add the 95% confidence band
                           risk.table = TRUE,      # Add risk table
                           tables.height = 0.2,
                           tables.theme = theme_cleantable(),
                           palette = "jco",
                           ggtheme = theme_bw(),
                           title = subset)
  print(p_surv_cut)
  csplots[[cg]] <-  p_surv_cut
  cg = cg + 1
}

length(csplots)
x1 = ceiling(sqrt(length(csplots)))
y1 = x1

all_plot <- arrange_ggsurvplots(csplots,
                                print = F,
                                ncol = x1, nrow = y1,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
# all_plot 
x2=5*x1
y2=5*y1
ggsave(all_plot, #path = prefix,
       filename = paste0(pro, 'all_cut_point_survival_plot.pdf'),width = x2,height = y2)





