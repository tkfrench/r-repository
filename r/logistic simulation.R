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
p=0.1; q=0.2; rho=0.45; n=1; n.sim <- 1000
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


#Hospital 1: Generate 1000 sample x,y cases, 0% leakage
#=====================================
set.seed(17);u <- a(rho, p, q);y <- floor((u-1)/2);x <- 1 - u %% 2
Hosp1 <- data.frame(hosp=1,x=x,y=y)
Hosp1$y_same_Hosp <-Hosp1$y * sample(c(rep(1,10)),n.sim,replace=T) # 0% leakage
Hosp1$Leakage = rep(mean(Hosp1$y) - mean(Hosp1$y_same_Hosp),1000)

#Hospital 2: Generate 1000 sample x,y cases, 10% leakage
#=====================================
set.seed(18);u <- a(rho, p, q);y <- floor((u-1)/2);x <- 1 - u %% 2
Hosp2 <- data.frame(hosp=2,x=x,y=y)
Hosp2$y_same_Hosp <-Hosp2$y * sample(c(rep(1,9),rep(0,1)),n.sim,replace=T) #10% leakage
Hosp2$Leakage = rep(mean(Hosp2$y) - mean(Hosp2$y_same_Hosp),1000)

#Hospital 3: Generate 1000 sample x,y cases, 20% leakage
#=====================================
set.seed(19);u <- a(rho, p, q);y <- floor((u-1)/2);x <- 1 - u %% 2
Hosp3 <- data.frame(hosp=3,x=x,y=y)
Hosp3$y_same_Hosp <-Hosp3$y * sample(c(rep(1,8),rep(0,2)),n.sim,replace=T) #20% leakage
Hosp3$Leakage = rep(mean(Hosp3$y) - mean(Hosp3$y_same_Hosp),1000)

#Hospital 4: Generate 1000 sample x,y cases, 40% leakage
#=====================================
set.seed(20);u <- a(rho, p, q);y <- floor((u-1)/2);x <- 1 - u %% 2
Hosp4 <- data.frame(hosp=4,x=x,y=y)
Hosp4$y_same_Hosp <-Hosp4$y * sample(c(rep(1,6),rep(0,4)),n.sim,replace=T) #40% leakage
Hosp4$Leakage = rep(mean(Hosp4$y) - mean(Hosp4$y_same_Hosp),1000)

# Model using all data
#=====================================
Hosp_All <- rbind(Hosp1, Hosp2, Hosp3, Hosp4)
 m1 <- glm(y ~ x , data = Hosp_All, family = binomial(link='logit'))
 pROC::roc(response = Hosp_All$y, predictor = fitted(m1))
 coef(summary(m1))[,'Pr(>|z|)']
 m1$coefficients
 b0=m1$coefficients[1]
 b1=m1$coefficients[2]


# fit all hospitals, evaluted using all captured outcomes.
#-----------------------------------------
 roc0 <- pROC::roc(response = Hosp_All$y,           predictor = fitted(m1))
 roc0
 Obs_rate = sum(Hosp_All$y)/nrow(Hosp_All)
 exp_rate = sum(exp(b0+b1*Hosp_All$x)/(1+exp(b0+b1*Hosp_All$x)))/nrow(Hosp_All)
 oe_ratio = Obs_rate/exp_rate 
 oe_ratio  # should = 1.00

#fit same model but evaluated assuming captured outcomes are same hospital only 
#-----------------------------------------
 roc1 <- pROC::roc(response = Hosp_All$y_same_Hosp, predictor = fitted(m1))
 roc1
 Obs_rate_same = sum(Hosp_All$y_same_Hosp)/nrow(Hosp_All)
 exp_rate = sum(exp(b0+b1*Hosp_All$x)/(1+exp(b0+b1*Hosp_All$x)))/nrow(Hosp_All)
 oe_ratio_same = Obs_rate_same/exp_rate 
 oe_ratio_same

# Results by hospital
#-----------------------------------------
Hosp_All %>% group_by(hosp) %>%
             summarize(Obs_rate=mean(y),
                       exp_rate=mean(exp(b0+b1*x)/(1+exp(b0+b1*x))),
                       OE_Ratio=Obs_rate/exp_rate,
                       Obs_rate_same=mean(y_same_Hosp),
                       OE_Ratio_same=Obs_rate_same/exp_rate)



 # Is there a difference in mean expected rate 
 # among event cases by location:same vs other facility 
 #-----------------------------------------

Hosp_All %>% group_by(hosp,y_same_Hosp) %>%
  filter(y==1) %>%
  summarize(exp_rate=mean(exp(b0+b1*x)/(1+exp(b0+b1*x))))          





