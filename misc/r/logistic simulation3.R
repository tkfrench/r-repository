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
set.seed(55)
hospitals       <- 51
hosp            <- c(0:(hospitals-1))
hospLeakageRate <- seq(0,0.50,0.50/(hospitals-1)) # Readmissions % at non-index facilities (0 to 50%)
hospEnctrs      <- rep(30,hospitals)            # Each hosp with with 1000 cases, change if desired
p               <- 0.1                            # Single risk factor 'x' prevalance     
q               <- 0.2                            # Overall readmission event rate 
rho             <- 0.45;                          # Spearman rank (rho) correlation between y and x 
n.sim           <- sum(hospEnctrs)
n               <- 1 # don't change
test<-list(hospLeakageRate =hospLeakageRate ,hosp=hosp)

for (i in 1:iterations) {
  u <- a(rho, p, q) # Randomize
  cohort <- data.frame(iteration=rep(i,n.sim),
                     hosp = rep(c(hosp),c(hospEnctrs)),
                     lkg  = rep(c(hospLeakageRate),c(hospEnctrs)),
                     x     = 1 - u %% 2,
                     y     = floor((u-1)/2)) 
  m1 <- glm(y ~ x , data = cohort, family = binomial(link='logit'))
  #pROC::roc(response = cohort$y, predictor = fitted(m1))
  #coef(summary(m1))[,'Pr(>|z|)']
  m1$coefficients
  b0=m1$coefficients[1]
  b1=m1$coefficients[2]
  cohort$expected <- exp(b0+b1*cohort$x)/(1+exp(b0+b1*cohort$x))
  # Initialize vectors
  if(i==1){  
    smry <- cohort %>% group_by(hosp, iteration) %>%
          summarize(Obs_rate=mean(y)
                 , Exp_rate =mean(expected)
                 , Obs_rate_same=mean(y*(1-lkg))
                 , OE_Ratio=mean(y)/mean(expected)
                 , OE_Ratio_same=mean(y*(1-lkg))/mean(expected)
                 , OE_Ratio_delta=mean(y)/mean(expected) - mean(y*(1-lkg))/mean(expected)
                 )
  }
  if(i>1){
    hold <- cohort %>% group_by(hosp, iteration) %>%
      summarize(Obs_rate=mean(y)
                , Exp_rate =mean(expected)
                , Obs_rate_same=mean(y*(1-lkg))
                , OE_Ratio=mean(y)/mean(expected)
                , OE_Ratio_same=mean(y*(1-lkg))/mean(expected)
                , OE_Ratio_delta=mean(y)/mean(expected) - mean(y*(1-lkg))/mean(expected)
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
        , OE_Ratio_same_avg =mean(OE_Ratio_same)
        , OE_Ratio_same_se  =sd(OE_Ratio_same)/sqrt(cases)
        , OE_Ratio_delta_avg=mean(OE_Ratio_delta)
        , OE_Ratio_delta_se  =sd(OE_Ratio_delta)/sqrt(cases)
          ) 

require(ggplot2)
ggplot(final, aes(x = hosp, y = OE_Ratio_same_avg)) +
  geom_point(size = 1)+
  geom_errorbar(aes(ymax = OE_Ratio_same_avg+2*OE_Ratio_same_se, ymin = OE_Ratio_same_avg-2*OE_Ratio_same_se))+
  geom_line(data = final,aes(x = hosp, y = OE_Ratio_same_avg), stat = "identity") +
  scale_y_continuous(limits = c(0.4, 1.05)) +
  xlab("Hosp (and leakage rate %)") +
  ylab('OE Ratio (Couning Other Admits Only')

require(ggplot2)
ggplot(final, aes(x = hosp, y = OE_Ratio_delta_avg)) +
  geom_point(size = 1)+
  geom_errorbar(aes(ymax = OE_Ratio_delta_avg+2*OE_Ratio_delta_se, ymin = OE_Ratio_delta_avg-2*OE_Ratio_delta_se))+
  geom_line(data = final,aes(x = hosp, y = OE_Ratio_delta_avg), stat = "identity") +
  scale_y_continuous(limits = c(0, 0.6)) +
  xlab("Hosp (and leakage rate %)") +
  ylab('OE Ratio (Couning Other Admits Only')



require(ggplot2)
ggplot(smry, aes(x = hosp, y = OE_Ratio_delta-hosp/100, group=hosp)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-.3, .3)) 



plot(smry$hosp, smry$OE_Ratio_same)
plot(smry$hosp, smry$OE_Ratio_delta)
plot(final$hosp, final$OE_Ratio_same_avg)
plot(final$hosp, final$OE_Ratio_delta_avg)
boxplot(OE_Ratio_same~hosp,data=smry, main="OE Ratios by varying Leakage rates", xlab="Hosp",ylab='OE Ratio')

require(ggplot2)
ggplot(smry, aes(x = hosp, y = OE_Ratio_same )) +
  geom_point(size = 1) +
  geom_line(data = final,aes(x = hosp, y = OE_Ratio_same_avg), stat = "identity")


require(ggplot2)
ggplot(smry, aes(x = hosp, y = OE_Ratio_delta, group=hosp)) +
  geom_boxplot() 






ggplot(final, aes(x = hosp, y = OE_Ratio_delta_avg)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymax = OE_Ratio_delta_avg+2*OE_Ratio_delta_se, ymin = OE_Ratio_delta_avg-2*OE_Ratio_delta_se)) +
  geom_line(data = final,aes(x = hosp, y = OE_Ratio_delta_avg), stat = "identity")
