rm(list=ls())
#library(dplyr)
library(plyr)
library(tidyverse)
library(qicharts2)
#------------------------------------------------------------------------
stratumVars         <- c("year","surgtype", "fracture")
BinaryOutcomeVars   <- c("AnyComplicationFlag", "Readmit30Flag")
IntervalOutcomeVars <- c("PTExpScore", "CompositeScore")
AllOutcomeVars      <- append(BinaryOutcomeVars, IntervalOutcomeVars)

# Generate sample cohort Assign:
#  - Entity type designations: region, ministry, provider
#  - Stratifiers: year, surgtype, fracture
#  - Outcomes Binary type(2, one rare); Interval type (2)
# --------------------------------------------------------------
caseNum <- 50000
provmnum =sample(1:1000,caseNum,replace=T); minNum = provmnum %%24+1; regNum = minNum %%5 +1
mydata <- data.frame(
  provider   =paste0("p",provmnum),     # calssify among 1000 providers
  ministry   =paste0("h",minNum),       # calssify among  24 Ministries
  region     =paste0("R",regNum),       # calssify among   5 Regions
  year       =sample(c(2016,2017,2018),caseNum,replace=T),
  mo         =sample(1:12,caseNum,replace=T),
  surgtype   =rbinom(caseNum,1, 0.40),   # Stratifier
  fracture   =rbinom(caseNum,1, 0.60),   # Stratifier
  AnyComplicationFlag =rbinom(caseNum,1, 0.01), # Binary outcome, Make it rare!
  AnyComplicationExp  =runif(caseNum,.001, .500), # Binary outcome, Make it rare!
  Readmit30Flag       =rbinom(caseNum,1, 0.30), # Binary outcome
  PTExpScore          =rnorm(caseNum,68,8),     # Interval outcome
  CompositeScore      =rnorm(caseNum,75,5)      # Interval toutcome
)
rm(caseNum, provmnum, minNum, regNum)
mydata$yr_mo <- 12*(mydata$year-min(mydata$year)) + mydata$mo
mydata$obs_exp <- mydata$AnyComplicationFlag - mydata$AnyComplicationExp

qicharts2::qic(y=obs_exp,  x= yr_mo, data  = mydata, chart = 'xbar')




# ------ fisher exact test ------------------------
# set muex = FALSE to compare vs System (cases aren't mutually exclusive)
#  T. French 4/20/2018
Fisher2by2.test <- function(num1,den1,num2,den2,muex=F){
  if(muex==FALSE){t<-matrix(c(num1, den1-num1, num2-num1, den2-num2-(den1-num1)),nrow=2)}
  else{           t<-matrix(c(num1, den1-num1, num2, den2-num2),nrow=2)}
  # test if too small for Chi Square
  if(min(t) <= 5 | sum(t)<=20){p.value = fisher.test(t)$p.value}
  else{p.value = chisq.test(t)$p.value}
  return(round(p.value,digits=5))}

# -Welch test -
# Assumes unequal variance  T. French 4/20/2018
welchsmryttest <- function(m1,m2,sd1,sd2,num1,num2,twosided=TRUE,alpha=0.05){
  se <- sqrt(sd1*sd1/num1+sd2*sd2/num2); t <- (m1-m2)/se
  df <-(sd1*sd1/num1+sd2*sd2/num2)^2 / (sd1^4/(num1^2*(num1-1))  + sd2^4/(num2^2*(num2-1)))
  if (twosided) {pvalue=2*(1-pt(abs(t),df=df))}
  else {pvalue=   1-pt(abs(t),df=df)}
  status <- ifelse(pvalue <= alpha ,ifelse((m1-m2) <0, "Lower","Higher"),"NS")
  #return <- list(pvalue=pvalue, t=t, df=df, status=status,alpha)
  return(round(pvalue,digits=5))}

# Funtion to generate list of all combinations of N variables [N choose k, k=1,..n]
#--------------------------------------------------------------------------------------
genVarCombos <-function(){listIn<-list()
for(i in 1:length(stratumVars)) {listIn[[length(listIn)+1]] <- combn(stratumVars,i)}
varCombos<-list()
for (j in 1:length(listIn)){
  for (k in 1:length(listIn[[j]][1,])){
    varCombos[[length(varCombos)+1]] <- list(c(unlist(listIn[[j]][,k])))}}
return(varCombos)}

varCombos<-genVarCombos()

#Summary results (mean, sd, n)
#     for all entities (Sys, rgn, min, prov),
#     for all outcomes
#     by all combinations of stratifiers (N chosen 1,2,...N at a time)
#----------------------------------------------------------------------------------------
systemStats<-data.frame();
regionStats<-data.frame();
ministryStats<-data.frame();
providerStats<-data.frame()
systemTest<-data.frame()
for(z in 1:length(varCombos)){
  the.list <- unlist(varCombos[(z)])
  temp  <-mydata %>% group_by_at(vars(the.list))%>%  summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA')
  systemStats<-rbind.fill(systemStats, temp )
  temp  <-mydata %>% group_by_at(vars(the.list,"region"))%>%  summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='region', entityID= region, ministry='NA',provider='NA')
  regionStats<-rbind.fill(regionStats, temp )
  temp  <-mydata %>% group_by_at(vars(the.list,"region","ministry"))%>% summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='ministry', entityID= ministry,  provider='NA')
  ministryStats<-rbind.fill(ministryStats, temp )
  temp  <-mydata %>% group_by_at(vars(the.list,"region","ministry","provider"))%>% summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='provider', entityID= provider)
  providerStats<-rbind.fill(providerStats, temp )

  temp  <-mydata %>% group_by_at(vars(the.list))%>%  summarise_at(vars(AllOutcomeVars),
   funs('mean_sys'=mean,  'sum_sys'=sum, 'sd_sys'=sd, 'n_sys'=sum(!is.na(.))), na.rm = TRUE)%>%
    mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA')
  systemTest<-rbind.fill(systemTest, temp )

  #systemTest summary results set for use with stat tests
  systemTest<-mydata %>%  group_by_at(vars(stratumVars))%>%
    summarise_at(vars(AllOutcomeVars), funs('mean_sys'=mean,  'sum_sys'=sum, 'sd_sys'=sd, 'n_sys'=sum(!is.na(.))), na.rm = TRUE)%>%
    mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA')
}
allStats<- rbind(systemStats, regionStats, ministryStats, providerStats)


#Entity vs system p-values for all [Binary outcomes], all entities, all stratum
#---------------------------------------------------------------------------
tic=proc.time()[3]
for(this.outcm in BinaryOutcomeVars){
  m  <-allStats   %>% select(stratumVars,contains(this.outcm),'region','ministry','provider')
  s  <-systemTest %>% select(stratumVars,contains(this.outcm)) %>%
    select(stratumVars,ends_with("_mean_sys"),ends_with("_sum_sys"),ends_with("_sd_sys"),ends_with("_n_sys"))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(Fisher2by2.test,",
                   "ms$",this.outcm,"_sum,","ms$",this.outcm,"_sum_sys,",
                   "ms$",this.outcm,"_n,",  "ms$",this.outcm,"_n_sys)" )
  eval(parse(text = eval_str))
  ms<- ms[, (colnames(ms) %in% c(stratumVars,'region','ministry','provider',paste0(this.outcm,"_pvalue") ))]
  allStats <- merge(allStats, ms, by = c(stratumVars,'region','ministry','provider')) # add p-values to aggregates
}
toc<-proc.time()[3] - tic; toc


#Entity vs system p-values for all [Interval outcomes], all entities, all stratum
#--------------------------------------------------------------------------------
tic=proc.time()[3]
for(this.outcm in IntervalOutcomeVars){
  m  <-allStats   %>% select(stratumVars,contains(this.outcm),'region','ministry','provider')
  s  <-systemTest %>% select(stratumVars,contains(this.outcm)) %>%
    select(stratumVars,ends_with("_mean_sys"),ends_with("_sum_sys"),ends_with("_sd_sys"),ends_with("_n_sys"))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(welchsmryttest,",
                   "ms$",this.outcm,"_mean,","ms$",this.outcm,"_mean_sys,",
                   "ms$",this.outcm,"_sd,",  "ms$",this.outcm,"_sd_sys,",
                   "ms$",this.outcm,"_n,",   "ms$",this.outcm,"_n_sys)" )
  eval(parse(text = eval_str))
  ms<- ms[, (colnames(ms) %in% c(stratumVars,'region','ministry','provider',paste0(this.outcm,"_pvalue") ))]
  allStats  <-merge(allStats , ms, by = c(stratumVars,'region','ministry','provider')) # add p-values to  Aggregates
}
toc<-proc.time()[3] - tic; toc

#Sort
allStats <- allStats[,order(names(allStats))]
allStats <- allStats[order(allStats$entityType, allStats$region, allStats$ministry, allStats$provider),]
#Reorder
allStats <- allStats %>% select(entityType, region, ministry, provider, entityID,  stratumVars, everything()) #reorder
head(allStats,5)

