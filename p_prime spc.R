rm(list=ls())
#library(dplyr)
#library(plyr)
library(tidyverse)
library(ggplot2)

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
  provider   =paste0("p",provmnum),     # calssify among 1000 providers
  ministry   =paste0("h",minNum),       # calssify among  24 Ministries
  region     =paste0("R",regNum),       # calssify among   5 Regions
  year       =sample(c(2016,2017,2018),caseNum,replace=T),
  month      =sample(1:12,caseNum,replace=T),
  surgtype   =rbinom(caseNum,1, 0.40),   # Stratifier
  fracture   =rbinom(caseNum,1, 0.60),   # Stratifier
  AnyComplicationFlag =rbinom(caseNum,1, 0.08), # Binary outcome, Make it rare!
  Readmit30Flag       =rbinom(caseNum,1, 0.30), # Binary outcome
  PTExpScore          =rnorm(caseNum,68,8),     # Interval outcome
  CompositeScore      =rnorm(caseNum,75,5)      # Interval toutcome
)
mydata$yrmo_ndx <- 12*(mydata$year-min(mydata$year)) + mydata$month


p_hat = mean(mydata$AnyComplicationFlag)
s_hat = sqrt(p_hat*(1-p_hat)/length(mydata$AnyComplicationFlag))
p_data <- mydata %>% group_by(yrmo_ndx) %>%
       summarize(num   = sum(AnyComplicationFlag),
                 denom = n(),
                 p   = num/denom,
                 s   = sqrt(p_hat*(1-p_hat)/denom),
                 cl  = p_hat,
                 lcl = max(0,p - 3*s),
                 ucl = min(1,p + 3*s),
                 
                 z_i   = (p-p_hat)/s
                # sigma_z = mean(abs(diff(z_i)), na.rm = TRUE) / 1.128,
                # s_prime = s*sigma_z,
                # lcl_prime = max(0,p - 3*s_prime),
               #  ucl_prime = max(0,p + 3*s_prime)
       )
sigma_z = mean(abs(diff(p_data$z_i)), na.rm = TRUE) / 1.128
s_prime = p_data$s*sigma_z
p_data$lcl_prime <- max(0,p_data$p - 3*s_prime)
p_data$ucl_prime <- max(0,p_data$p + 3*s_prime)

ggplot(data = p_data, aes(x=yrmo_ndx)) + 
  geom_line(aes(y = p), colour="blue") + 
  geom_line(aes(y = cl), colour="blue") + 
  geom_line(aes(y = ucl), colour="grey") + 
  geom_line(aes(y = lcl), colour="grey")   

ggplot(data = p_data, aes(x=yrmo_ndx)) + 
  geom_line(aes(y = p ), colour="blue") + 
  geom_line(aes(y = cl), colour="blue") + 
  geom_line(aes(y = ucl_prime), colour="grey") + 
  geom_line(aes(y = lcl_prime), colour="grey")  