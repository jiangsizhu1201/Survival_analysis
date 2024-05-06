rm(list=ls())
library(survival)
library(survminer) 
library(ggstatsplot) 
proj='tcga-luad'
load(file = paste0(proj,".for_survival.rdata") )

# 1. prepare data for coxph----
## 批量生存分析使用coxph回归方法
## http://www.sthda.com/english/wiki/cox-proportional-hazasc_dataset-model
exprSet[1:4,1:4]
phe=meta
head(phe)
mySurv <- with(phe, Surv(time, event))
survival_dat=phe

# 测试一些基因看看是否有生存分析意义
gene= as.numeric(exprSet['PRC1',])
gene= as.numeric(exprSet['DMD',])
survival_dat$group=ifelse(gene>median(gene),'high','low')  
fit <- survfit(Surv(time, event) ~ group,
               data = survival_dat) 
ggsurvplot(fit,data = survival_dat, 
                   pval = T, #在图上添加log rank检验的p值
                   # pval.method = TRUE,#添加p值的检验方法
                   risk.table = TRUE, 
                   risk.table.y.text = F,
                   xlab = "Time in years", #x轴标题
                   # xlim = c(0, 10), #展示x轴的范围
                   # break.time.by = 1, #x轴间隔
                   size = 1.5, #线条大小
                   ggtheme = theme_ggstatsplot(),
                   palette="nejm" )


cox_results <-apply(exprSet , 1 , function(gene){
  # gene= as.numeric(exprSet[1,]) 
  group=ifelse(gene>median(gene),'high','low') 
  group=factor(  group,levels = c('low','high'))
  if( length(table(group))<2)
    return(NULL)
  survival_dat <- data.frame(group=group,# stage=phe$stage,
                             stringsAsFactors = F)
  
  m=coxph(mySurv ~ group, 
          # mySurv ~  stage+ group,  # 如果有交叉变量
          data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, 
                     p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouphigh',])
  
})

# 2. specify the value----
cox_results=t(cox_results)
head(cox_results)
cox_results=cox_results[order(cox_results[,'HR'],decreasing = T),]

table(cox_results[,4]<0.01)
table(cox_results[,4]<0.05)
cox_results['DMD',]
save(cox_results, 
     file = 'batch_cox_results.Rdata')
load(     file = 'batch_cox_results.Rdata')


