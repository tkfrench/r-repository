rm(list=ls())
#library(dplyr)
#library(plyr)
library(tidyverse)

#------------------------------------------------------------------------
stratumVars         <- c("year","surgtype", "fracture")
BinaryOutcomeVars   <- c("AnyComplicationFlag", "Readmit30Flag")
IntervalOutcomeVars <- c("PTExpScore", "CompositeScore")
AllOutcomeVars      <- append(BinaryOutcomeVars, IntervalOutcomeVars)

# Generate sample cohort; assigning:
#  - Entity type: region, ministry, provider
#  - Stratifiers: year, surgtype, fracture
#  - Outcomes Binary type(2, one rare); Interval type (2)
# --------------------------------------------------------------
caseNum <- 5000
provmnum =sample(1:1000,caseNum,replace=T); minNum = provmnum %%24+1; regNum = minNum %%5 +1
mydata <- data.frame(
  cohort     =sample(c(replicate(20, "A"), replicate(5, "B"), "C"), caseNum, replace=T),
  provider   =paste0("p",provmnum),     # calssify among 1000 providers
  ministry   =paste0("h",minNum),       # calssify among  24 Ministries
  region     =paste0("R",regNum),       # calssify among   5 Regions
  year       =sample(c(2016,2017,2018),caseNum,replace=T),
  month      =sample(1:12,caseNum,replace=T),
  surgtype   =rbinom(caseNum,1, 0.40),   # Stratifier
  fracture   =rbinom(caseNum,1, 0.60),   # Stratifier
  outlier    =rbinom(caseNum,1, 0.50),   # Stratifier
  AnyComplicationFlag =rbinom(caseNum,1, 0.01), # Binary outcome, Make it rare!
  Readmit30Flag       =rbinom(caseNum,1, 0.30), # Binary outcome
  PTExpScore          =rnorm(caseNum, 0.68, 0.08),     # Interval outcome
  LTRxpScore          =rnorm(caseNum, 0.78, 0.10)     # Interval outcome
  #CompositeScore      =rnorm(caseNum,075,5)      # Interval toutcome
)
rm(caseNum, provmnum, minNum, regNum)

mydata$comp <- (2*mydata$PTExpScore + 3*mydata$LTRxpScore)/(2+3)

# creatt PTExpScore2 that includes missing values 
mydata$PTExpScore2 <- ifelse(mydata$outlier==1,NA,mydata$PTExpScore)

by_cohort <- mydata %>% group_by(cohort) %>% 
  summarize(Sys_cases            = n(),
            Sys_mean = mean(PTExpScore2,       na.rm = TRUE),
            Sys_var  = var(PTExpScore2,        na.rm = TRUE),
            Sys_n    = sum(!is.na(PTExpScore2),na.rm = TRUE),
            LTRxpScore_mean = mean(LTRxpScore,       na.rm = TRUE),
            LTRxpScore_var  = var(LTRxpScore,       na.rm = TRUE),
            comp_mean = mean(comp,       na.rm = TRUE),
            comp_var = var(comp,       na.rm = TRUE))

by_cohort$comp_mean2 = (2*by_cohort$Sys_mean+ 3*by_cohort$LTRxpScore_mean)/(2+3)
by_cohort$comp_var2  = (4*by_cohort$Sys_var*(by_cohort$Sys_n-1) + 9*by_cohort$LTRxpScore_var*(by_cohort$Sys_n-1))/
                        (4+9+ by_cohort$Sys_n-1+by_cohort$Sys_n-1)

by_cohort_ministry <- mydata  %>% group_by(cohort, ministry) %>% 
  summarize(cases            = n(),
            mean = mean(PTExpScore2,       na.rm = TRUE),
            var  = var(PTExpScore2,        na.rm = TRUE),
            n    = sum(!is.na(PTExpScore2),na.rm = TRUE)) 

test <- merge(by_cohort,by_cohort_ministry, all=TRUE )

test$partial_var2 <- test$n * test$Sys_var2 / test$Sys_n
test
test$varRatio <- test$Sys_var/test$var
test$varRatio <- ifelse(test$varRatio>15000 | test$n < 5 , NA, test$varRatio)

ggplot(test, aes(varRatio, n)) +geom_point()
ggplot(test, aes(x=varRatio)) +geom_histogram()
ggplot(test, aes(x=varRatio)) +geom_density()
summary(test$varRatio)

mydata %>%
  mutate(region  = factor(region),
         ministry= factor(ministry),
         provider= factor(provider)) %>%
  map(~ mean(.x)) %>%
  map_dfr(~ broom::tidy(.), .id = 'source') %>%
 #mutate(p.value = round(p.value, 5))
  mutate(p.value = 0.0)



#by_provider <- mydata %>%
#  group_by(region, ministry, provider ) %>%
#  nest() %>%
#  map(data, mean()) %>%


# ------ 2 by 2 test of pre-aggregated data  ----------------------
# set muex = FALSE to compare vs System (cases aren't mutually exclusive); T. French 4/20/2018
smry2By2.test <- function(num1,den1,num2,den2,muex=F){
  if(muex==FALSE){t<-matrix(c(num1, den1-num1, num2-num1, den2-num2-(den1-num1)),nrow=2)}
  else{           t<-matrix(c(num1, den1-num1, num2, den2-num2),nrow=2)}
  # test if too small for Chi Square
  if(min(t) <= 5 | sum(t)<=20){p.value = fisher.test(t)$p.value}
  else{p.value = chisq.test(t)$p.value}
  return(round(p.value,digits=5))}

# ------ Welch test for pre-aggregated data  ----------------------
# Assumes unequal variance  T. French 4/20/2018
smryWelcht.test <- function(m1,m2,sd1,sd2,num1,num2,twosided=TRUE,alpha=0.05){
  se <- sqrt(sd1*sd1/num1+sd2*sd2/num2); t <- (m1-m2)/se
  df <-(sd1*sd1/num1+sd2*sd2/num2)^2 / (sd1^4/(num1^2*(num1-1))  + sd2^4/(num2^2*(num2-1)))
  if (twosided) {pvalue=2*(1-pt(abs(t),df=df))}
  else {pvalue=   1-pt(abs(t),df=df)}
  status <- ifelse(pvalue <= alpha ,ifelse((m1-m2) <0, "Lower","Higher"),"NS")
  #return <- list(pvalue=pvalue, t=t, df=df, status=status,alpha)
  return(round(pvalue,digits=5))}


# Generate list of all combinations of N stratum variables [N choose K, K=1,..N]; usning an 'anonymous function'
  stratCombos <-(function(x){listIn<-list()
  for(i in 1:length(x)) {listIn[[length(listIn)+1]] <- combn(x,i)}
  combos<-list()
  for (j in 1:length(listIn)){
    for (k in 1:length(listIn[[j]][1,])){
    combos[[length(combos)+1]] <- list(c(unlist(listIn[[j]][,k])))}}
  return(combos)})(stratumVars)

#Summary results (mean, sum, sd, n)
#     for all entities (Sys, rgn, min, prov),
#     for all outcomes
#     by all combinations of stratifiers (N chosen 1,2,...N at a time)
#----------------------------------------------------------------------------------------
  tic=proc.time()[3]
systemStats<-data.frame(); regionStats<-data.frame(); ministryStats<-data.frame(); providerStats<-data.frame()
systemTest <-data.frame()
for(z in 1:length(stratCombos)){
  the.list <- unlist(stratCombos[(z)])

  systemStats<-rbind.fill(systemStats,
    mydata %>% group_by_at(vars(the.list))%>%  summarise_at(vars(AllOutcomeVars),
    funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA'))

  regionStats<-rbind.fill(regionStats,
    mydata %>% group_by_at(vars(the.list,"region"))%>%  summarise_at(vars(AllOutcomeVars),
    funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='region', entityID= region, ministry='NA',provider='NA'))

  ministryStats<-rbind.fill(ministryStats,
    mydata %>% group_by_at(vars(the.list,"region","ministry"))%>% summarise_at(vars(AllOutcomeVars),
    funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='ministry', entityID= ministry,  provider='NA'))

  providerStats<-rbind.fill(providerStats,
    mydata %>% group_by_at(vars(the.list,"region","ministry","provider"))%>% summarise_at(vars(AllOutcomeVars),
    funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
    mutate(entityType  ='provider', entityID= provider))

  #systemTest summary results set for use with stat tests, includes system 25, 50, 75th percentiles
  systemTest<-rbind.fill(systemTest,
    mydata %>% group_by_at(vars(the.list))%>% summarise_at(vars(AllOutcomeVars),
    funs('mean_sys'=mean,  'sum_sys'=sum, 'sd_sys'=sd, 'n_sys'=sum(!is.na(.)),
    q25=quantile(., probs=0.25), q75=quantile(., probs=0.75), q50=quantile(., probs=0.50)), na.rm = TRUE)%>%
    mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA'))
}
allStats<- rbind(systemStats, regionStats, ministryStats, providerStats)
 toc<-proc.time()[3] - tic; toc

#Replace "NA" with 'All' for stratum with this.strat levels (these are the don't care elements)
for(this.strat in stratumVars){
  eval_str=paste0("allStats$",this.strat," <- as.character(allStats$",this.strat,")");	eval(parse(text = eval_str))
  eval_str=paste0("allStats$",this.strat,"[is.na(allStats$",this.strat,")]<-'All'");  	eval(parse(text = eval_str))

  eval_str=paste0("systemTest$",this.strat," <- as.character(systemTest$",this.strat,")"); eval(parse(text = eval_str))
  eval_str=paste0("systemTest$",this.strat,"[is.na(systemTest$",this.strat,")]<-'All'");  	eval(parse(text = eval_str))
}

#Entity vs system p-values for all [Binary outcomes], all entities, all stratum combinations
#---------------------------------------------------------------------------
tic=proc.time()[3]
for(this.outcm in BinaryOutcomeVars){
  m  <-allStats   %>% select(stratumVars,contains(this.outcm),'region','ministry','provider')
  s  <-systemTest %>% select(stratumVars,contains(this.outcm)) %>%
    select(stratumVars,ends_with("_mean_sys"),ends_with("_sum_sys"),ends_with("_sd_sys"),ends_with("_n_sys"))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(smry2By2.test,",
                   "ms$",this.outcm,"_sum,","ms$",this.outcm,"_sum_sys,",
                   "ms$",this.outcm,"_n,",  "ms$",this.outcm,"_n_sys)" )
  eval(parse(text = eval_str))
  ms<- ms[, (colnames(ms) %in% c(stratumVars,'region','ministry','provider',paste0(this.outcm,"_pvalue") ))]
  allStats <- merge(allStats, ms, by = c(stratumVars,'region','ministry','provider')) # add p-values to aggregates
}
toc<-proc.time()[3] - tic; toc

#Entity vs system p-values for all [Interval outcomes], all entities, all stratum combinations
#--------------------------------------------------------------------------------
tic=proc.time()[3]
for(this.outcm in IntervalOutcomeVars){
  m  <-allStats   %>% select(stratumVars,contains(this.outcm),'region','ministry','provider')
  s  <-systemTest %>% select(stratumVars,contains(this.outcm)) %>%
    select(stratumVars,ends_with("_mean_sys"),ends_with("_sum_sys"),ends_with("_sd_sys"),ends_with("_n_sys"))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(smryWelcht.test,",
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

allStats %>%
  filter(entityType=='provider') %>%
  select(provider, surgtype, fracture, CompositeScore_mean,CompositeScore_pvalue, CompositeScore_n) %>%
  mutate(surgtype_fracture=paste0(surgtype,"_",fracture)) %>%
  mutate(signif =cut(CompositeScore_pvalue, breaks=c(0.0,0.01,0.05), include.lowest=TRUE)) %>%
  #spread(term, estimate) %>%
  ggplot(aes(CompositeScore_mean, surgtype_fracture)) +
  geom_jitter(aes(colour= signif, size=CompositeScore_n),alpha=0.3,height=.25)

allStats %>%
  filter(entityType=='ministry') %>%
  select(ministry, surgtype, fracture, CompositeScore_mean,CompositeScore_pvalue, CompositeScore_n) %>%
  mutate(surgtype_fracture=paste0(surgtype,"_",fracture)) %>%
  mutate(signif =cut(CompositeScore_pvalue, breaks=c(0.0,0.01,0.05,1), include.lowest=TRUE)) %>%
  #spread(term, estimate) %>%
  ggplot(aes(CompositeScore_mean, surgtype_fracture)) +
  geom_jitter(aes(colour= signif, size=CompositeScore_n),alpha=0.3,height=.25)


allStats %>%
  filter(entityType=='region') %>%
  select(region, surgtype, fracture, CompositeScore_mean,CompositeScore_pvalue, CompositeScore_n) %>%
  mutate(surgtype_fracture=paste0(surgtype,"_",fracture)) %>%
  mutate(signif =cut(CompositeScore_pvalue, breaks=c(0.0,0.01,0.05,1))) %>%
  #spread(term, estimate) %>%
  ggplot(aes(CompositeScore_mean, surgtype_fracture)) +
  geom_jitter(aes(colour= signif, size=CompositeScore_n),alpha=0.3,height=.25)

mu <- 100
s <- 50
n <- 5
nsim <- 10000 # number of simulations
# theoretical standard error
s / sqrt(n)
# simulation of experiment and the standard deviations of their means
y <- replicate( nsim, mean( rnorm(n, mu, s) ) )
sd(y)
