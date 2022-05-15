
rm(list=ls())
library(tidyverse) # only needed for data simulation

#======= begin Simulate Raw data =====================
n = 10000
df <- data.frame(ndx = seq(1:n), facility = round(1+51*rbeta(n, 2, 2)), 
                 dschrg_dt = sample(seq(as.Date('2018/01/01'), as.Date('2019/12/31'), by="day"), n, replace=T),
                 #post = sample(0:1, n, replace=T),
                 mort_exp = rbeta(n, 1, 5), rdmt_exp = rbeta(n, 1, 5), cost = rnorm(n,mean=5, sd=2))
  df$mort_obs  <- rbinom(n, size=1, df$mort_exp); df$rdmt_obs <- rbinom(n, size=1, df$rdmt_exp)
  df$cost_sqrd <- df$cost * df$cost
  df$post      <- ifelse(df$dschrg_dt > '2018/12/31', 1,0)
 # randomly inject missing values into x% of cases
  df$mort_obs[sample(seq(1:n), .01*n, replace=F)] <- rep(NA,.01*n); df$mort_exp[sample(seq(1:n), .01*n, replace=F)] <- rep(NA,.01*n)
  df$rdmt_obs[sample(seq(1:n), .03*n, replace=F)] <- rep(NA,.03*n); df$rdmt_exp[sample(seq(1:n), .03*n, replace=F)] <- rep(NA,.03*n)
  df$cost[sample(seq(1:n), .02*n, replace=F)]     <- rep(NA,.02*n); df$cost_sqrd[sample(seq(1:n),.02*n, replace=F)] <- rep(NA,.02*n)
df <- data.frame(ndx=df$ndx, facility=df$facility, dschrg_dt= df$dschrg_dt, post=df$post, 
                 ent_mort_obs = df$mort_obs, ent_mort_pre_obs = df$mort_obs*(1-df$post),ent_mort_post_obs = df$mort_obs*(df$post),
                 ent_mort_exp = df$mort_exp, ent_mort_pre_exp = df$mort_exp*(1-df$post),ent_mort_post_exp = df$mort_exp*(df$post),
                 ent_rdmt_obs = df$rdmt_obs, ent_rdmt_pre_obs = df$rdmt_obs*(1-df$post),ent_rdmt_post_obs = df$rdmt_obs*(df$post),
                 ent_rdmt_exp = df$rdmt_exp, ent_rdmt_pre_exp = df$rdmt_exp*(1-df$post),ent_rdmt_post_exp = df$rdmt_exp*(df$post),
                 ent_cost     = df$cost,     ent_cost_pre     = df$cost*(1-df$post),    ent_cost_post     = df$cost*(df$post),
                 ent_cost_sqrd =df$cost,     ent_cost_sqrd_pre= df$cost_sqrd*(1-df$post),    ent_cost_sqrd_post     = df$cost_sqrd*(df$post)) 
#======= end Simulate RAW data =====================
#======= begin Simulate AGGREGATE data =====================
aggByFac <- df %>% group_by(facility) %>%
              summarise_at(., 
                           vars(ent_mort_obs, ent_mort_pre_obs, ent_mort_post_obs, 
                                ent_mort_exp, ent_mort_pre_exp, ent_mort_post_exp,
                                ent_rdmt_obs, ent_rdmt_pre_obs, ent_rdmt_post_obs, 
                                ent_rdmt_exp, ent_rdmt_pre_exp, ent_rdmt_post_exp,
                                ent_cost,     ent_cost_pre,     ent_cost_post,
                                ent_cost_sqrd,ent_cost_sqrd_pre,ent_cost_sqrd_post),
                           funs(sum, n=sum(!is.na(.)) ), na.rm=T)
aggByFacYrMo <- df %>% mutate(month = format(dschrg_dt, "%m"), year = format(dschrg_dt, "%Y")) %>%
  group_by(facility, month, year) %>%
  summarise_at(., 
               vars(ent_mort_obs, ent_mort_pre_obs, ent_mort_post_obs, 
                    ent_mort_exp, ent_mort_pre_exp, ent_mort_post_exp,
                    ent_rdmt_obs, ent_rdmt_pre_obs, ent_rdmt_post_obs, 
                    ent_rdmt_exp, ent_rdmt_pre_exp, ent_rdmt_post_exp,
                    ent_cost,     ent_cost_pre,     ent_cost_post,
                    ent_cost_sqrd,ent_cost_sqrd_pre,ent_cost_sqrd_post),
               funs(sum, n=sum(!is.na(.)) ), na.rm=T)
# Clean up
rm(list=setdiff(ls(), c("df", "aggByFac", "aggByFacYrMo")))
#======= End Simulate aggByFacREGATE data =====================


## ************** ASYMPTOTIC TEST FOR THE RATIO OF TWO POISSON RATES, O/E'S, OR SIR's  ***************
ratiotest <- function( obs1, obs2, exp1, exp2, alpha=0.05){
  obs1  <- round(obs1); obs2 <- round(obs2)
  OE1   <- ifelse((obs1 >= 0 & exp1 > 0), obs1/exp1, NaN)
  OE2   <- ifelse((obs2 >= 0 & exp2 > 0), obs2/exp2, NaN)
  zstat <- ifelse((exp1>0 & exp2>0 & is.finite(OE1) & is.finite(OE2)),
                  2*(sqrt(obs1 + 0.375) - sqrt((exp1/exp2)*(obs2 + 0.375))) / sqrt(1 + (exp1/exp2)),NaN)
  p.value         <- ifelse((exp1>0 & exp2>0 & is.finite(OE1) & is.finite(OE2)),2*(1 - pnorm(abs(zstat))),NaN)
  signif_status   <- ifelse(p.value <= alpha & OE2 > OE1, "hi",
                            ifelse(p.value <= alpha & OE1 > OE2, "lo",'ns'))    
  p.value.str     <- paste0(substr(formatC(p.value, width = 6, digits=4, format='f'), 0,6 ))
  event_delta     <- round(obs2 - (OE1 * exp2))
  event_delta.str <- paste0(substr(formatC(event_delta, width = 6, digits=0, format='f'), 0,6 ))
  return(list(p.value     = p.value,     event_delta     = event_delta, 
              p.value.str = p.value.str, event_delta.str = event_delta.str,
              signif_status = signif_status))
} 
## ****************************************************************************************************

#  Example OE Ratio tests: Each entity vs rest
#---------------------------------------------------------
OE_v_rest <- function( obs, exp, alpha = 0.05){
  sys_obs_sum = sum(obs, rm.na=T); sys_exp_sum = sum(exp, rm.na=T)
  df <- data.frame(ent_obs_sum = obs, ent_exp_sum = exp, alpha = alpha,
                   rest_obs_sum = sys_obs_sum - obs,
                   rest_exp_sum = sys_exp_sum - exp)
  pre_v_rest_test <- ratiotest(df$rest_obs_sum, df$ent_obs_sum, df$rest_exp_sum, df$ent_exp_sum, alpha=alpha)
  paste0("|vrest.p|", pre_v_rest_test$p.value.str, "|vrestevent_delta|", pre_v_rest_test$event_delta.str,
         "|signif_status|", pre_v_rest_test$signif_status)
}

.arg1 = aggByFac$ent_mort_obs_sum 
.arg2 = aggByFac$ent_mort_exp_sum 
.arg3 = 0.05
OE_v_rest(.arg1, .arg2, .arg3)
#pre_v_rest_test <- OE_v_rest(.arg1, .arg2,.arg3)
#paste0("|vrest.p|", pre_v_rest_test$p.value.str, "|vrestevent_delta|", pre_v_rest_test$event_delta.str)

#sys_obs_sum = sum(.arg1, rm.na=T); sys_exp_sum = sum(.arg2, rm.na=T)
#df <- data.frame(ent_obs_sum = .arg1, ent_exp_sum = .arg2, alpha=.arg3,
#                 rest_obs_sum = sys_obs_sum - .arg1,
#                 rest_exp_sum = sys_exp_sum - .arg2)
#
#pre_v_rest_test <- ratiotest(df$rest_obs_sum, df$ent_obs_sum, df$rest_exp_sum, df$ent_exp_sum, alpha=alpha)
#paste0("|vrest.p|", pre_v_rest_test$p.value.str, "|vrestevent_delta|", pre_v_rest_test$event_delta.str)
# returns  strings that looks like:  "|vrest.p|0.5527|vrestevent_delta|    -3"

#-  Example OE Ratio tests: Each entity pre vs. post
#---------------------------------------------------------
.arg1 = aggByFac$ent_mort_pre_obs_sum 
.arg2 = aggByFac$ent_mort_post_obs_sum
.arg3 = aggByFac$ent_mort_pre_exp_sum 
.arg4 = aggByFac$ent_mort_post_exp_sum 
.arg5 = 0.05

#
df <- data.frame(ent_pre_obs_sum = .arg1, ent_post_obs_sum = .arg2,
                 ent_pre_exp_sum = .arg3, ent_post_exp_sum = .arg4,
                 alpha = .arg5)

ent_post_test   <- ratiotest(df$ent_pre_obs_sum, df$ent_post_obs_sum, df$ent_pre_exp_sum, df$ent_post_exp_sum, alpha=alpha)
paste0("|vpost.p|", ent_post_test$p.value.str, "|vpost_event_delta|", ent_post_test$event_delta.str)
# returns strings that looks like: "|vpost.p|0.1550|vpost_event_delta|    -6" 

#---- Example OE Ratio tests: BOTH entity vs rest & pre vs post
.arg1 = aggByFac$ent_mort_obs_sum 
.arg2 = aggByFac$ent_mort_pre_obs_sum 
.arg3 = aggByFac$ent_mort_post_obs_sum
.arg4 = aggByFac$ent_mort_exp_sum 
.arg5 = aggByFac$ent_mort_pre_exp_sum 
.arg6 = aggByFac$ent_mort_post_exp_sum
.arg7 = 0.05

#
sys_obs_sum = sum(.arg1, rm.na=T); sys_exp_sum = sum(.arg4, rm.na=T)
df <- data.frame(ent_obs_sum = .arg1, ent_pre_obs_sum = .arg2, ent_post_obs_sum = .arg3, 
                 ent_exp_sum = .arg4, ent_pre_exp_sum = .arg5, ent_post_exp_sum = .arg6,
                 rest_obs_sum = sys_obs_sum - .arg1,
                 rest_exp_sum = sys_exp_sum - .arg4,
                 alpha=.arg7)

pre_v_rest_test <- ratiotest(df$rest_obs_sum, df$ent_obs_sum, df$rest_exp_sum, df$ent_exp_sum, alpha=alpha)
ent_post_test   <- ratiotest(df$ent_pre_obs_sum, df$ent_post_obs_sum, df$ent_pre_exp_sum, df$ent_post_exp_sum, alpha=alpha)
paste0("|vrest.p|", pre_v_rest_test$p.value.str, "|vrest_event_delta|", pre_v_rest_test$event_delta.str,
       "|vpost.p|", ent_post_test$p.value.str,   "|vpost_event_delta|", ent_post_test$event_delta.str)
# returns strings that looks like:  "|vrest.p|0.5368|vrest_event_delta|     2|vpost.p|0.6397|vpost_event_delta|     2"

#==============================
WilsonCI <- function(num, den, alpha=0.05){
  estimate = num/den
  z = qnorm(1-alpha/2, 0, 1)
  A = (2*den*estimate + z^2)/(2*den + 2*z^2)
  B = (z*sqrt(z^2 + 4*den*estimate*(1 - estimate)))/(2*den + 2*z^2)
  output<-list(max(c(0,(A-B))), min(c(1, (A+B))))
  return(output)
}      #use mapply(WilsonCI, intake$num, intake$den, intake$alpha)

alpha = 0.05
Wt_mort = .5
Wt_rdmt = .5


## ************** Composite scoring through the Markov Chain Ranking Method  *************** 
dim = nrow(smry)
# ---- Mortality OE
mort_mat <- matrix(0, dim , dim) 
for (i in 1: dim){
  for (j in 1:dim){
    num1 = smry$ent_mort_obs[i]; num2 = smry$ent_mort_obs[j]
    den1 = smry$ent_mort_exp[i]; den2 = smry$ent_mort_exp[j]
    p    = ratiotest(num1, num2, den1, den2)$p.value
    mort_mat[i,j] = ifelse(i == j, 0.0, 
                             #  ifelse( p > alpha, 0.5,
                             #      ifelse( p < alpha & (num1/den1 > num2/den2), 1, 0)))                            
                           ifelse( (num1/den1 >  num2/den2), 0.5 + p/2, 0.5 - (1 - p)/2))
  }
}
# ---- Readmission OE
rdmt_mat <- matrix(0, dim , dim) 
for (i in 1: dim){
  for (j in 1:dim){
    num1 = smry$ent_rdmt_obs[i];  num2 = smry$ent_rdmt_obs[j]
    den1 = smry$ent_rdmt_exp[i];  den2 = smry$ent_rdmt_exp[j]
    p    = ratiotest(num1, num2, den1, den2)$p.value
    rdmt_mat[i,j] = ifelse(i == j, 0.0, 
                          # ifelse( p > alpha, 0.5,
                          #        ifelse( p < alpha & (num1/den1 > num2/den2), 1, 0)))                            
                        ifelse( (num1/den1 >  num2/den2), 0.5 + p/2, 0.5 - (1 - p)/2))
  }
}

# Row normalize each outcome NxN matrix
mort_mat_norm <- mort_mat * matrix( rep( 1/rowSums(mort_mat, na.rm=T)), dim, dim)
rdmt_mat_norm <- rdmt_mat * matrix( rep( 1/rowSums(rdmt_mat, na.rm=T)), dim, dim)
# Replace outcome NxN matrix row entries with 1/n if rowsum == 0
mort_mat_norm[rowSums( mort_mat_norm != 0, na.rm=T) == 0, ] <- rep( 1/dim, dim)
rdmt_mat_norm[rowSums( rdmt_mat_norm != 0, na.rm=T) == 0, ] <- rep( 1/dim, dim)
# Combine all (weighted) outcome NxN matricies as matrix V
V <- (mort_mat_norm * Wt_mort) + (rdmt_mat_norm * Wt_rdmt) 
# Return eigenvector corresponding to dominate eigenvalue
e <- eigen(t(V), symmetric = FALSE)
smry$MarkovComposite <- abs(Re(e$vectors[,1]))

plot(smry$ent_rdmt_obs/smry$ent_rdmt_exp, smry$MarkovComposite)
cor(smry$ent_rdmt_obs/smry$ent_rdmt_exp, smry$MarkovComposite)

plot(smry$ent_mort_obs/smry$ent_mort_exp, smry$MarkovComposite)
cor(smry$ent_mort_obs/smry$ent_mort_exp, smry$MarkovComposite)

plot(smry$ent_rdmt_obs/smry$ent_rdmt_exp+smry$ent_mort_obs/smry$ent_mort_exp, smry$MarkovComposite)
cor(smry$ent_rdmt_obs/smry$ent_rdmt_exp+smry$ent_mort_obs/smry$ent_mort_exp, smry$MarkovComposite)

# ----------------------------------------------------------------------
#  WELCH T TEST   
# ----------------------------------------------------------------------
# Function: p.value comparing paired averages w unequal variances
welch_t.pvalue <- function( m1, m2, sd1, sd2, num1, num2, alpha = 0.05, twosided = TRUE){
  se <- ifelse((num1 > 0 & num2 > 0), sqrt(sd1*sd1/num1+sd2*sd2/num2), NaN)
  t  <- ifelse((is.finite(m1) & is.finite(m2) & is.finite(se)), (m1-m2)/se, 0.0)
  df <- ifelse((num1 > 0 & num2 > 0), (sd1*sd1/num1+sd2*sd2/num2)^2 / (sd1^4/(num1^2*(num1-1))  + sd2^4/(num2^2*(num2-1))), NaN)
  if (twosided) {p.value = ifelse((is.finite(t) & is.finite(df)), 2*(1-pt(abs(t),df=df)), 1.0)} 
  else {p.value = ifelse((is.finite(t) & is.finite(df)), 1-pt(abs(t),df=df), 1.0)}  
  return(p.value)
}
