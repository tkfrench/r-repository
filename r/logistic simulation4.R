##################################################################
# Routine for generating sample random dataset for testing Logistic model
#  Generates two variables - y: response variable; x: risk factor
#  q: ~Prevalence of y (set by user)
#  p: ~Prevalence of x (set by user)
#  rho: Spearman rank correlation
#       used as strength of association of y with x, (set by user).
#       Higher the abs(rho) => higher ROC for glm(y ~ x,family = binomial(link='logit'))
# Four hospital sets of 1,000 cases generated with varying leakage rates% (0,10,20,40)
##################################################################
rm(list=ls())
library(tidyverse)
library(reshape2)
library(pROC)

# This snippet generates the data, called once for each hospital
# All hospitals: y prevalence (q) ~ 20%; x prevalence (p) ~10%, N=1000
#-----------------------------------------

# Compute Pr(0,0) from rho, p=Pr(X=1), and q=Pr(Y=1).
a <- function(rho, p, q) {
   a.0  <- rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
   prob <- c(`(0,0)`=a.0, `(1,0)`=1-q-a.0, `(0,1)`=1-p-a.0, `(1,1)`=a.0+p+q-1)
   if (min(prob) < 0) {
     print(prob)
     stop("Error: a probability is negative.")
   }
   u <- sample.int(4, n.sim * n, replace=TRUE, prob=prob)
   return(u)
}
#-----------------------------------------
iterations = 100
staticEnctern = 100
set.seed(21)
hospitals       <- 51
hospid          <- c(0:(hospitals-1))
#hospLeakageRate <- seq(0,0.50,0.50/(hospitals-1)) # Readmissions % at non-index facilities (0 to 50%)
hospEnctrs      <- rep(staticEnctern,hospitals)            # Each hosp with with 1000 cases, change if desired
p               <- 0.1                            # Single risk factor 'x' prevalance     
q               <- 0.3                            # Overall readmission event rate 
rho             <- 0.45;                          # Spearman rank (rho) correlation between y and x 
n.sim           <- sum(hospEnctrs)
n               <- 1 # don't change


for (i in 1:iterations) {
  u <- a(rho, p, q) # Randomize
  cohort <- data.frame(iteration=rep(i,n.sim),
                     hosp = rep(c(hospid),c(hospEnctrs)),
                     #lkg  = rep(c(hospLeakageRate),c(hospEnctrs)),
                     lkg  = rep(c(hospid/100),c(hospEnctrs)),
                     x     = 1 - u %% 2,
                     y     = floor((u-1)/2),
                     y_ndxhosp = floor((u-1)/2)*unlist(lapply(rep(c(1-hospid/100),c(hospEnctrs)), function(x) rbinom(n=1, size=1,prob=x)))
                     )
  m1 <- glm(y ~ x , data = cohort, family = binomial(link='logit'))
  #pROC::roc(response = cohort$y, predictor = fitted(m1))
  #coef(summary(m1))[,'Pr(>|z|)']
  m1$coefficients
  b0=m1$coefficients[1]
  b1=m1$coefficients[2]
  cohort$expected <- exp(b0+b1*cohort$x)/(1+exp(b0+b1*cohort$x))
  cohort$AUC <-0*cohort$y + pROC::roc(response = cohort$y, predictor = fitted(m1))$auc[1]
  # Initialize vectors
  if(i==1){  
    smry <- cohort %>% group_by(hosp, iteration) %>%
          summarize(Obs_rate=mean(y)
                 , Exp_rate      = mean(expected)
                 , Obs_rate_same = mean(y_ndxhosp)
                 , OE_Ratio      = mean(y)/mean(expected)
                 , OE_Ratio_same = mean(y_ndxhosp )/mean(expected)
                 , OE_Ratio_same_adj       = Obs_rate_same/ mean((1-lkg)*expected)
                 , OE_Ratio_same_adj_delta = mean(y)/mean(expected) - Obs_rate_same/ mean((1-lkg)*expected)
                 , Rate_Ratio    = (Obs_rate_same/ mean((1-lkg)*expected)) / (mean(y)/mean(expected))
                 , AUC    = mean(AUC)
                 )
  }
  if(i>1){
    hold <- cohort %>% group_by(hosp, iteration) %>%
      summarize(Obs_rate=mean(y)
                , Exp_rate       = mean(expected)
                , Obs_rate_same  = mean(y_ndxhosp)
                , OE_Ratio       = mean(y)/mean(expected)
                , OE_Ratio_same  = mean(y_ndxhosp )/mean(expected)
                , OE_Ratio_same_adj       = Obs_rate_same/ mean((1-lkg)*expected)
                , OE_Ratio_same_adj_delta = mean(y)/mean(expected) - Obs_rate_same/ mean((1-lkg)*expected)
                , Rate_Ratio    = (Obs_rate_same/ mean((1-lkg)*expected)) / (mean(y)/mean(expected))
                , AUC    = mean(AUC)
      )
    smry<-rbind(smry,hold)    
  } 

}
  
final<-smry %>% group_by(hosp) %>%
summarize(cases = n()
        , Obs_rate_avg      =mean(Obs_rate)
        , Exp_rate_avg      =mean(Exp_rate)
        , Obs_rate_same_avg =mean(Obs_rate_same)
        , OE_Ratio_avg      =mean(OE_Ratio)
        , OE_Ratio_25       =quantile(OE_Ratio, .27)
        , OE_Ratio_50       =quantile(OE_Ratio, .50)
        , OE_Ratio_75       =quantile(OE_Ratio, .75)
        , OE_Ratio_same_avg =mean(OE_Ratio_same)
        , OE_Ratio_same_50  =quantile(OE_Ratio_same, .50)
        , OE_Ratio_same_se  =sd(OE_Ratio_same)/sqrt(cases)
        , OE_Ratio_same_adj_avg =mean(OE_Ratio_same_adj)
        , OE_Ratio_same_adj_05 = quantile(OE_Ratio_same_adj, .05)
        , OE_Ratio_same_adj_25 = quantile(OE_Ratio_same_adj, .25) 
        , OE_Ratio_same_adj_50 = quantile(OE_Ratio_same_adj, .50) 
        , OE_Ratio_same_adj_75 = quantile(OE_Ratio_same_adj, .75)        
        , OE_Ratio_same_adj_95 = quantile(OE_Ratio_same_adj, .95)
        , Rate_Ratio_05 = quantile(Rate_Ratio, .05)
        , Rate_Ratio_25 = quantile(Rate_Ratio, .25) 
        , Rate_Ratio_50 = quantile(Rate_Ratio, .50) 
        , Rate_Ratio_75 = quantile(Rate_Ratio, .75)        
        , Rate_Ratio_95 = quantile(Rate_Ratio, .95)  
          ) 


dd_sub = final[,c("hosp","OE_Ratio_50","OE_Ratio_same_adj_50","OE_Ratio_same_50")]
dd = melt(dd_sub, id=c("hosp"))
ggplot(dd) + geom_line(aes(x=hosp, y=value, colour=variable)) +
  xlab("Hosp (leakage rate%)") + ylab("OE Ratios") +
  scale_colour_manual(values=c("green","blue","red")) +
  labs(title=paste0("Median OE Ratios: 51 Ministries with 0% to 50% Leakage Rates"
                    ,"\nRandomized with:"
                    ,"\n  Readmision Prevalence ~ ",round(100*q,0),"%"
                    ,"\n  Single Dependent Variable Prevalence ~ ",round(100*p,0),"%"
                    ,"\n  Logistic Model: Rho(y,x)=0.45 => Model ROC ~0.65"
                    ,"\n\n   Cases per Hosp: ",staticEnctern,",  Iterations ",iterations))


dd_sub = final[,c("hosp","Rate_Ratio_25","Rate_Ratio_50","Rate_Ratio_75")]
dd = melt(dd_sub, id=c("hosp")) 
ggplot(dd) + geom_line(aes(x=hosp, y=value, colour=variable))+ 
  scale_colour_manual(values=c("blue","green","red")) + 
  xlab("Hosp (leakage rate%)") + ylab("Rate Ratio") +
  geom_line(aes(x=hosp, y=value, colour=variable)) 

ggplot(data=smry,aes(y=AUC)) +
geom_boxplot()

ggplot(data=smry,aes(y=AUC,x=iteration )) +
  geom_jitter(width = 0.5, size=1,alpha=.10)

AUCs <- smry %>% group_by(iteration)%>%
                 summarize(AUC=mean(AUC))
ggplot(data=AUCs,aes(x=AUC )) +
  geom_density()
ggplot(data=AUCs,aes(x=AUC )) +
  geom_histogram()
