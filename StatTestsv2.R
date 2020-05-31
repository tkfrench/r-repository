rm(list=ls())
require(dplyr)
library(plyr)
# -----  Welch test function-----------------------
welchsmryttest <- function(m1,m2,sd1,sd2,num1,num2,twosided=TRUE,alpha=0.05){
  se <- sqrt(sd1*sd1/num1+sd2*sd2/num2); t <- (m1-m2)/se
  df <-(sd1*sd1/num1+sd2*sd2/num2)^2 / (sd1^4/(num1^2*(num1-1))  + sd2^4/(num2^2*(num2-1)))
  if (twosided) {pvalue=2*(1-pt(abs(t),df=df))} 
  else {pvalue=   1-pt(abs(t),df=df)}  
  return(round(pvalue,digits=5))
}#-------------------------------------------------
# ------ fisher exact test ------------------------
#   set muex = FALSE to compare vs system (cases aren't mutually exclusive)
Fisher2by2.test <- function(num1,den1,num2,den2,muex=F){
  if(muex==FALSE){t<-matrix(c(num1, den1-num1, num2-num1, den2-num2-(den1-num1)),nrow=2)} 
  else{           t<-matrix(c(num1, den1-num1, num2, den2-num2),nrow=2)}
  # test if too big for Fisher test
  if(((num1+num2)*((den1-num1)+(den2-num2)))< 2^31 - 1){p.value = fisher.test(t)$p.value}
  else{p.value = chisq.test(t)$p.value}
  return(round(p.value,digits=5))
} #--------------------------------------------------
# ------  Standarized Rate Ratio test ( for comparing OE ratios, SIRs, etc.) ----------
SMRtest <- function(x1,x2,t1,t2){
  # assumes Indirectly Standardized Rates (ISR) follow a Poisson ditribution
  if(t1>0 & t2>0){p.value <- poisson.test(c(x1,x2),c(t1,t2))$p.value}
  else{p.value=NaN}
  return(p.value)
  #SMRtest(8, 5, 3, 3.1 )
}#-------------------------------------------------

# Generate 100,000 random encounters 
# with [5 regions],[24 ministires], [2 startifiers] & [4 outcomes]

stratumVars         <- c("gender",   "age_gt65")
BinaryOutcomeVars   <- c("Readmit", "AnyComplication", "Mortality")
IntervalOutcomeVars <- c("LOS", "PtExp")
AllOutcomeVars      <- append(BinaryOutcomeVars, IntervalOutcomeVars)
CaseCount=100000
rndmnum =sample(1:100,CaseCount,replace=T)
mydata<-data.frame( 
  region   = paste0("R",1+ rndmnum %%5),       # calssify among   5 regions  
  ministry = paste0("h",1+ rndmnum %%24),      # calssify among  24 Ministries
  provider = paste0("p",1+ rndmnum %%100),     # calssify among 100 providers    
  gender   = sample(c("M", "F"),CaseCount,replace=T), # Add statifier
  age_gt65 = sample(c("Y", "N"),CaseCount,replace=T), # Add statifier
  LOS      = c(rnorm(n=CaseCount, mean=50, sd=10)),   # Add Interval outcome
  PtExp    = c(rnorm(n=CaseCount, mean=60, sd=9)),    # Add Interval outcome
  Readmit  = sample(c(0,1),CaseCount,replace=T),      # Add Binary outcome
  Mortality= sample(c(0,1),CaseCount,replace=T),      # Add Binary outcome
  AnyComplication= sample(c(0,1),CaseCount,replace=T))# Add Binary outcome

#Entity type summary results (mean, sd, n) by stratifies for numeric outcomes 1 thru N
systemStats  <-mydata %>% group_by_at(vars(stratumVars))%>%           summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>% 
               mutate(entityType  ='system', entityID= 'system', region='NA',ministry='NA',provider='NA')
regionStats  <-mydata %>% group_by_at(vars(stratumVars,"region"))%>%  summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
               mutate(entityType  ='region', entityID= region, region=region,ministry='NA',provider='NA')
ministryStats<-mydata %>% group_by_at(vars(stratumVars,"region","ministry"))%>%summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
               mutate(entityType  ='ministry', entityID= region, region=region,ministry=ministry,provider='NA')
providerStats<-mydata %>% group_by_at(vars(stratumVars,"region","ministry","provider"))%>%summarise_at(vars(AllOutcomeVars), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE) %>%
               mutate(entityType  ='provider', entityID= region, region=region,ministry=ministry,provider=provider)
allStats<- rbind(systemStats,regionStats,ministryStats,providerStats)


#ministry v system p-values for [Interval outcomes]
for(this.outcm in IntervalOutcomeVars){
#for(this.outcm in c('outcome1')){  
  this.outcm.list    =c(paste0(this.outcm,"_mean"),     paste0( this.outcm,"_sd"),     paste0(this.outcm,"_n"))
  this.outcm.list_sys=c(paste0(this.outcm,"_mean_sys"), paste0( this.outcm,"_sd_sys"), paste0(this.outcm,"_n_sys"))
  m<-data.frame(subset(allStats,select=c(stratumVars,'region','ministry','provider',     this.outcm.list)))
  s<-data.frame(subset(systemStats,  select=c(stratumVars,this.outcm.list)))
  names(s) <- sub("_mean", "_mean_sys", names(s))
  names(s) <- sub("_sd",   "_sd_sys",   names(s))
  names(s) <- sub("_n",    "_n_sys",    names(s))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(welchsmryttest,",
                   "ms$",this.outcm.list[1],",","ms$",this.outcm.list_sys[1],",",  # means
                   "ms$",this.outcm.list[2],",","ms$",this.outcm.list_sys[2],",",    # SDs
                   "ms$",this.outcm.list[3],",","ms$",this.outcm.list_sys[3],")" )   # Ns
  eval(parse(text = eval_str)) 
  ms<- ms[, (colnames(ms) %in% c(stratumVars,'region','ministry','provider',paste0(this.outcm,"_pvalue") ))]
  allStats  <-merge(allStats , ms, by = c(stratumVars,'region','ministry','provider')) # add p-values to  Aggregates
}


#Entity v system p-values for [Binary outcomes]
for(this.outcm in BinaryOutcomeVars){
  this.outcm.list    =c(paste0(this.outcm,"_sum"),     paste0( this.outcm,"_n"))
  this.outcm.list_sys=c(paste0(this.outcm,"_sum_sys"), paste0( this.outcm,"_n_sys"))
  m<-data.frame(subset(allStats,select=c(stratumVars,'region','ministry','provider',     this.outcm.list)))
  s<-data.frame(subset(systemStats,  select=c(stratumVars,this.outcm.list)))
  names(s) <- sub("_sum", "_sum_sys", names(s))
  names(s) <- sub("_n",    "_n_sys",  names(s))
  ms <-merge(m, s, by = c(stratumVars))
  eval_str<-paste0("ms$",this.outcm,"_pvalue= mapply(Fisher2by2.test,",
                   "ms$",this.outcm.list[1],",","ms$",this.outcm.list_sys[1],",",  # num
                   "ms$",this.outcm.list[2],",","ms$",this.outcm.list_sys[2],")" ) # den
  eval(parse(text = eval_str)) 
  ms<- ms[, (colnames(ms) %in% c(stratumVars,'region','ministry','provider',paste0(this.outcm,"_pvalue") ))]
  allStats  <-merge(allStats , ms, by = c(stratumVars,'region','ministry','provider')) # add p-values to  Aggregates
}  
rm(list=setdiff(ls(), "allStats"))
allStats <- allStats %>% select(entityType, region, ministry, provider, entityID, everything()) #reorder
allStats <-allStats[order(allStats$entityType, allStats$region, allStats$ministry, allStats$provider),]
head(allStats,5)

