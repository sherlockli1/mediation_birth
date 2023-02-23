##################################################################################################################
# Code used to conduct analyses for "Contributions of the social and physical environments to racial and ethnic  #
# disparities in birth outcomes in Los Angeles County from 2017 to 2019                                          #
# Analytic approach developed by: Sherlock Li                                                                    #
# Code by: Sherlock Li                                                                                           #
# Last Update: 2/23/2023                                                                                         #
##################################################################################################################

require(haven)
library("dplyr")
library(tidyr)
library(gtsummary)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(tidyverse)
library(caret)
library(leaps)
library(lme4)
library(foreach)
library(doParallel)
library(hrbrthemes)
library(mice)
library(FactoMineR)
library(glmm)

final_birth_1<-read.csv("finaldata.csv")

#Subset into White vs Black
final_birth_bw<-final_birth_1[which(final_birth_1$PGB_MULTIRACE_116==1|final_birth_1$PGB_MULTIRACE_116==2),] #change here
final_birth_bw$m_race<-ifelse(final_birth_bw$PGB_MULTIRACE_116==2,1,0)

#Subset into White vs Hispanic
final_birth_his<-final_birth_1[which(final_birth_1$PGB_MULTIRACE_116==1|final_birth_1$PGB_MULTIRACE_116==8),] #change here
final_birth_his$m_race<-ifelse(final_birth_his$PGB_MULTIRACE_116==8,1,0)

#Correlation Matrix
myvars<-c("below_highschool","lingisolated","belowpoverty","unemployment",
          "NDVI","water","block_pm")
final_birth_cor<-final_birth_1[,myvars]
COR_mediator<-cor(final_birth_cor,use = "complete.obs",method="pearson")
write.csv(COR_mediator,"cor_mediator.csv")

#Table 1
final_birth_1 %>%
  dplyr::select(PGB_MULTIRACE_116,
         ptb,tlbw,
         m_edu_cat,payment_prenatal,m_birthplace_US,
         CHILD_SEX_14,smoking,parity,m_age_cat,BMI_cat,index_prenatal,marital)%>%
  tbl_summary(
    by = PGB_MULTIRACE_116,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    missing_text = "(Missing)"
  )

final_birth_1 %>%
  dplyr::select(PGB_MULTIRACE_116,
                below_highschool,lingisolated,belowpoverty,unemployment,
                water,block_pm,NDVI)%>%
  tbl_summary(
    by = PGB_MULTIRACE_116,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    missing_text = "(Missing)"
  )

#Black vs White
#PTB
data1<-final_birth_bw[,c("m_race",
                         "ptb",
                         "n_ses_score",
                         "env_score_3",
                         "m_edu_cat",
                         "payment_prenatal",
                         "m_birthplace_US",
                         "CHILD_SEX_14" ,
                         "smoking",
                         "parity",
                         "m_age_cat",
                         "BMI_cat",
                         "index_prenatal",
                         "marital",
                         "GISJOIN")]
#TLBW
data2<-final_birth_bw[,c("m_race",
                         "tlbw",
                         "n_ses_score",
                         "env_score_3",
                         "m_edu_cat",
                         "payment_prenatal",
                         "m_birthplace_US",
                         "CHILD_SEX_14" ,
                         "smoking",
                         "parity",
                         "m_age_cat",
                         "BMI_cat",
                         "index_prenatal",
                         "marital",
                         "GISJOIN")]


#Hispanic vs White
#PTB
data1a<-final_birth_his[,c("m_race",
                           "ptb",
                           "n_ses_score",
                           "env_score_3",
                           "m_edu_cat",
                           "payment_prenatal",
                           "m_birthplace_US",
                           "CHILD_SEX_14" ,
                           "smoking",
                           "parity",
                           "m_age_cat",
                           "BMI_cat",
                           "index_prenatal",
                           "marital",
                           "GISJOIN")]

#TLBW
data2a<-final_birth_his[,c("m_race",
                           "tlbw",
                           "n_ses_score",
                           "env_score_3",
                           "m_edu_cat",
                           "payment_prenatal",
                           "m_birthplace_US",
                           "CHILD_SEX_14" ,
                           "smoking",
                           "parity",
                           "m_age_cat",
                           "BMI_cat",
                           "index_prenatal",
                           "marital",
                           "GISJOIN")]


process_data=function(data){
  colnames(data) <- c('a',
                      'y',
                      'm1a',
                      'm1b',
                      'c1',
                      'c2',
                      'c3',
                      'c4',
                      'c5',
                      'c6',
                      'c7',
                      'c8',
                      'c9',
                      'c10',
                      'ID')
  data <- data[complete.cases(data), ]
  data$c1<-as.factor(data$c1)
  data$c2<-as.factor(data$c2)
  data$c3<-as.factor(data$c3)
  data$c4<-as.factor(data$c4)
  data$c5<-as.factor(data$c5)
  data$c6<-as.factor(data$c6)
  data$c7<-as.factor(data$c7)
  data$c8<-as.factor(data$c8)
  data$c9<-as.factor(data$c9)
  data$c10<-as.factor(data$c10)
  
  return(data)}

data1<-process_data(data1)
data2<-process_data(data2)

data1a<-process_data(data1a)
data2a<-process_data(data2a)


#Parallel Mediators
pse_para=function(data){
  data=as.data.frame(data)
  m1a_model =lm(m1a ~ a,data=data)
  m1b_model =lm(m1b ~ a,data=data)
  y_model   =lmer(y ~ a+m1a+m1b+c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+a:m1a+a:m1b+(1|ID),data=data)
  c1d=data$c1
  c2d=data$c2
  c3d=data$c3
  c4d=data$c4
  c5d=data$c5
  c6d=data$c6
  c7d=data$c7
  c8d=data$c8
  c9d=data$c9
  c10d=data$c10
  IDd=data$ID
  
  #G computation for Path Specific Effects
  #m1(x=1)
  m_1a_a=predict(m1a_model,data.frame(a=1))
  m_1b_a=predict(m1b_model,data.frame(a=1))
  
  #m1(x=0)
  m_1a_b=predict(m1a_model,data.frame(a=0))
  m_1b_b=predict(m1b_model,data.frame(a=0))

  #Pure Indirect
  #Vary Mediation 1a fix exposure to 1
  ma=predict(y_model,data.frame(
    m1a=m_1a_a,
    m1b=m_1b_b,
    a=1,
    c1=c1d,
    c2=c2d,
    c3=c3d,
    c4=c4d,
    c5=c5d,
    c6=c6d,
    c7=c7d,
    c8=c8d,
    c9=c9d,
    c10=c10d,
    ID=IDd))
  #Vary 1b
  mb=predict(y_model,data.frame(
    m1a=m_1a_b,
    m1b=m_1b_a,
    a=1,
    c1=c1d,
    c2=c2d,
    c3=c3d,
    c4=c4d,
    c5=c5d,
    c6=c6d,
    c7=c7d,
    c8=c8d,
    c9=c9d,
    c10=c10d,
    ID=IDd))
  #reference for indirect effect
  mc=predict(y_model,data.frame(
    m1a=m_1a_b,
    m1b=m_1b_b,
    a=1,
    c1=c1d,
    c2=c2d,
    c3=c3d,
    c4=c4d,
    c5=c5d,
    c6=c6d,
    c7=c7d,
    c8=c8d,
    c9=c9d,
    c10=c10d,
    ID=IDd))
  md=predict(y_model,data.frame(
    m1a=m_1a_b,
    m1b=m_1b_b,
    a=0,
    c1=c1d,
    c2=c2d,
    c3=c3d,
    c4=c4d,
    c5=c5d,
    c6=c6d,
    c7=c7d,
    c8=c8d,
    c9=c9d,
    c10=c10d,
    ID=IDd))
  #Total effect
  me=predict(y_model,data.frame(
    m1a=m_1a_a,
    m1b=m_1b_a,
    a=1,
    c1=c1d,
    c2=c2d,
    c3=c3d,
    c4=c4d,
    c5=c5d,
    c6=c6d,
    c7=c7d,
    c8=c8d,
    c9=c9d,
    c10=c10d,
    ID=IDd))
  #Total Indirect
  tie1a_md=mean(ma)-mean(mc)
  tie1b_md=mean(mb)-mean(mc)
  #Total Effect
  total_md=mean(me)-mean(md)
  return(c(tie1a_md,tie1b_md,total_md))
}

############################################
bootstrap=function(data){
  r=sample(1:n1,n1,replace=TRUE)
  ndata=data[r,]
  return(ndata)
}

#Process Parallel Mediator Results
process_tie_result=function(result){
  result<-as.data.frame(t(result))
  colnames(result)<-c("tie1a_md","tie1b_md","total_md")
  result_mean_1a<-round(mean(result$tie1a_md),digits = 4)
  result_mean_1b<-round(mean(result$tie1b_md),digits = 4)
  result_mean_tot<-round(mean(result$total_md),digits = 4)
  
  result_ci_1a<-round(quantile(result$tie1a_md,probs=c(0.025,0.975)),digits = 4)
  result_ci_1b<-round(quantile(result$tie1b_md,probs=c(0.025,0.975)),digits = 4)
  result_ci_tot<-round(quantile(result$total_md,probs=c(0.025,0.975)),digits = 4)
  
  result_final_mean<-cbind(result_mean_1a,
                           result_mean_1b,
                           result_mean_tot)
  
  result_final_ci<-cbind(result_ci_1a,
                         result_ci_1b,
                         result_ci_tot)
  
  result_final<-rbind(result_final_mean,result_final_ci)
  result_final<-as.data.frame(t(result_final))
  result_final$change<-round((result_final[,1])/(result_final[3,1])*100,digits = 2)
  result_final$change_L<-round((result_final[,2])/(result_final[3,1])*100,digits = 2)
  result_final$change_U<-round((result_final[,3])/(result_final[3,1])*100,digits = 2)
  return(result_final)}


set.seed(seed=1234)

#n1=50000 or 100000
n2=200
result=array(data = NA,dim=c(3,n2))


n1=nrow(data1)
for (i in 1:n2) {result[,i]=pse_para(bootstrap(data1))}
result1<-as.data.frame(result)
result1<-process_tie_result(result1)

n1=nrow(data2)
for (i in 1:n2) {result[,i]=pse_para(bootstrap(data2))}
result2<-as.data.frame(result)
result2<-process_tie_result(result2)

n1=nrow(data1a)
for (i in 1:n2) {result[,i]=pse_para(bootstrap(data1a))}
result1a<-as.data.frame(result)
result1a<-process_tie_result(result1a)

n1=nrow(data2a)
for (i in 1:n2) {result[,i]=pse_para(bootstrap(data2a))}
result2a<-as.data.frame(result)
result2a<-process_tie_result(result2a)


result_all_bw<-cbind(result1,
                     result2)
result_all_his<-cbind(result1a,
                     result2a)

write.csv(result_all_bw,"results/result_all_bw_newmodel3.csv")
write.csv(result_all_his,"results/result_all_his_newmodel3.csv")


















