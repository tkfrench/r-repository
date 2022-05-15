rm(list=ls())
#require(dplyr)
library(tidyverse)
library(broom)

#------------------------------------------------------------------------
# Generate sample cohort Assign:
#  - Entity type designations: region, ministry, provider
#  - Stratifiers: year, surgtype, fracture
#  - Outcomes Binary type(2, one rare); Interval type (2)
# --------------------------------------------------------------
stratumVars         <- c("year","surgtype", "fracture")
BinaryOutcomeVars   <- c("AnyComplicationFlag", "Readmit30Flag")
IntervalOutcomeVars <- c("PTExpScore", "CompositeScore")
AllOutcomeVars      <- append(BinaryOutcomeVars, IntervalOutcomeVars)
caseNum <- 50000
 
provmnum =sample(1:1000,caseNum,replace=T); minNum = provmnum %%24+1; regNum = minNum %%5 +1
mydata <- tibble(stringsAsFactors=FALSE,
  provider   =paste0("p",provmnum),     # calssify among 1000 providers    
  ministry   =paste0("h",minNum),       # calssify among  24 Ministries
  region     =paste0("R",regNum),       # calssify among   5 Regions
  QtrSeq    =sample(c(1,2,3,4,5,6,7,8),caseNum,replace=T),  
  surgtype   =rbinom(caseNum,1, 0.40),   # Stratifier
  fracture   =rbinom(caseNum,1, 0.60),   # Stratifier
  AnyComplicationFlag =rbinom(caseNum,1, 0.01), # Binary outcome, Make it rare!
  AnyComplicationExp  =runif(caseNum,.001, 0.02), # Binary outcome, Make it rare!
  Readmit30Flag       =rbinom(caseNum,1, 0.30), # Binary outcome
  PTExpScore          =rnorm(caseNum,68,8),     # Interval outcome
  CompositeScore      =rnorm(caseNum,75,5)      # Interval toutcome
)
rm(caseNum, provmnum, minNum, regNum)

# What entity type to aggregate and test?
#---------------------------------
entity.type='ministry'
#entity.type='region'
#entity.type='provider'


if(entity.type=='ministry') {
    mydata$the.entity  = mydata$ministry
  } else if (entity.type=='region') {
    mydata$the.entity  = mydata$region
  } else  {
    mydata$the.entity  = mydata$provider
  }

#mydata$the.entity  = mydata$ministry
mydata$the.period  = mydata$QtrSeq

#======================================================
# Test of trend by Entity, Interval type, encounter level data:
mydata$the.outcome =  mydata$CompositeScore
#======================================================
#mydata$the.entity = mydata$ministry
dd= mydata %>% group_by(the.entity) %>%
    do({model = lm(the.outcome ~QtrSeq, data=.) # create model
        data.frame(tidy(model),                    # get coefficient info
        broom:::glance(model))})                   # get model info
dd %>% filter(term=='QtrSeq') %>% select(the.entity, estimate, p.value) 

#======================================================
# Test of trend by Entity, Binomial type, encounter level data - logit:
mydata$the.outcome =  mydata$Readmit30Flag
#======================================================
mydata$the.entity = mydata$ministry
dd= mydata %>% group_by(the.entity) %>%
  do({model = glm(the.outcome ~QtrSeq, data=., family=binomial(link = "logit")) # create model
  data.frame(tidy(model),                               # get coefficient info
             broom:::glance(model))})                   # get model info
dd %>% filter(term=='QtrSeq') %>% select(the.entity, estimate, p.value) 


#======================================================
mydata$the.outcome =  mydata$Readmit30Flag
# Test of trend by Entity, Binomial type, Summary level data - logit:
#======================================================
mydata$the.entity = mydata$ministry
PropByEntByPeriod <- mydata %>% group_by(the.entity, the.period) %>% 
  summarise(num =sum(the.outcome),
            denom=sum(!is.na(the.outcome)),
            prop =sum(the.outcome)/sum(!is.na(the.outcome))) 

dd= PropByEntByPeriod %>% group_by(the.entity) %>%
  do({model = glm(cbind(num,(denom-num)) ~ the.period, family=binomial(logit), data=.)
  data.frame(tidy(model),                    # get coefficient info
             broom:::glance(model))})                   # get model info
dd %>% filter(term=='the.period') %>% select(the.entity, estimate, p.value) 

#======================================================
# Test of trend by Entity, Binomial type, Summary level data - Poisson:
mydata$the.outcome =  mydata$Readmit30Flag
#======================================================
mydata$the.entity = mydata$ministry
PropByEntByPeriod <- mydata %>% group_by(the.entity, the.period) %>% 
  summarise(num =sum(the.outcome),
            denom=sum(!is.na(the.outcome)),
            prop =sum(the.outcome)/sum(!is.na(the.outcome))) 

dd= PropByEntByPeriod %>% group_by(the.entity) %>%
  do({model = glm(num ~ denom + the.period, family="poisson", data=.) # create model
  data.frame(tidy(model),                    # get coefficient info
             broom:::glance(model))})                   # get model info
dd %>% filter(term=='the.period') %>% select(the.entity, estimate, p.value) 

