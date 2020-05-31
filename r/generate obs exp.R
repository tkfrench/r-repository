
rm(list=ls())
library(tidyverse)
alpha = 0.05
Wt_mort = .5
Wt_rdmt = .5
#======= begin Simulate data =====================
n = 5000
df <- data.frame(ndx = seq(1:n), facility = sample(1:52, n, replace=T), post = sample(0:1, n, replace=T),
                 mort_exp = rbeta(n, 1, 5), rdmt_exp = rbeta(n, 1, 5))

df$mort_obs <- rbinom(n, size=1, df$mort_exp)
df$rdmt_obs <- rbinom(n, size=1, df$rdmt_exp)

sys_mort_obs <- sum(df$mort_obs); sys_mort_exp=sum(df$mort_exp)
sys_rdmt_obs <- sum(df$rdmt_obs); sys_rdmt_exp=sum(df$rdmt_exp)
# -- rollup aggregate results into smry
smry <- df %>% 
              group_by(facility) %>%
               summarise(ent_mort_obs = sum(mort_obs), ent_mort_exp = sum(mort_exp),
                         ent_rdmt_obs = sum(rdmt_obs), ent_rdmt_exp = sum(rdmt_exp),
                         ent_mort_pre_obs  = sum((mort_obs*(1-post))), ent_mort_pre_exp  = sum(mort_exp*(1-post)),
                         ent_mort_post_obs = sum((mort_obs*(  post))), ent_mort_post_exp = sum(mort_exp*(  post)),
                         ent_rdmt_pre_obs  = sum((rdmt_obs*(1-post))), ent_rdmt_pre_exp  = sum(rdmt_exp*(1-post)),
                         ent_rdmt_post_obs = sum((rdmt_obs*(  post))), ent_rdmt_post_exp = sum(rdmt_exp*(  post)))  %>%
               mutate(rest_mort_obs = sys_mort_obs - ent_mort_obs,
                      rest_mort_exp = sys_mort_exp - ent_mort_exp,
                      rest_rdmt_obs = sys_rdmt_obs - ent_rdmt_obs,
                      rest_rdmt_exp = sys_rdmt_exp - ent_rdmt_exp)
#======= end Simulate data =====================

WilsonCI <- function(num, den, alpha=0.05){
  estimate = num/den
  z = qnorm(1-alpha/2, 0, 1)
  A = (2*den*estimate + z^2)/(2*den + 2*z^2)
  B = (z*sqrt(z^2 + 4*den*estimate*(1 - estimate)))/(2*den + 2*z^2)
  output<-list(max(c(0,(A-B))), min(c(1, (A+B))))
  return(output)
}      #use mapply(WilsonCI, intake$num, intake$den, intake$alpha)

## ************** ASYMPTOTIC TEST FOR THE RATIO OF TWO POISSON RATES, O/E'S, OR SIR's  ***************
ratiotest <- function( obs1, obs2, exp1, exp2, alpha=0.05){
  obs1  <- round(obs1); obs2 <- round(obs2)
  OE1   <- ifelse((obs1 >= 0 & exp1 > 0), obs1/exp1, NaN)
  OE2   <- ifelse((obs2 >= 0 & exp2 > 0), obs2/exp2, NaN)
  zstat <- ifelse((exp1>0 & exp2>0 & is.finite(OE1) & is.finite(OE2)),
                  2*(sqrt(obs1 + 0.375) - sqrt((exp1/exp2)*(obs2 + 0.375))) / sqrt(1 + (exp1/exp2)),NaN)
  p.value         <- ifelse((exp1>0 & exp2>0 & is.finite(OE1) & is.finite(OE2)),2*(1 - pnorm(abs(zstat))),NaN)
  p.value.str     <- paste0(substr(formatC(p.value, width = 6, digits=4, format='f'), 0,6 ))
  event_delta     <- round(obs2 - (OE1 * exp2))
  event_delta.str <- paste0(substr(formatC(event_delta, width = 6, digits=0, format='f'), 0,6 ))
  return(list(p.value     = p.value,     event_delta     = event_delta, 
              p.value.str = p.value.str, event_delta.str = event_delta.str))
} 
## ****************************************************************************************************
smry$OutputStr  <- rep('',nrow(smry))

pre_post_test   <- ratiotest(smry$rest_mort_obs, smry$ent_mort_obs, smry$rest_mort_exp, smry$ent_mort_exp, alpha=0.05)
ent_v_rest_test <- ratiotest(smry$ent_mort_pre_obs, smry$ent_mort_post_obs, smry$ent_mort_pre_exp, smry$ent_mort_post_exp, alpha=0.05)

smry$OutputStr  <- paste0(smry$OutputStr,'|pre_post_p.value|', pre_post_test$p.value.str)
smry$OutputStr  <- paste0(smry$OutputStr,'|v_rest_p.value|',   ent_v_rest_test$p.value.str)

smry$OutputStr  <- paste0(smry$OutputStr,'|pre_post_delta|',  pre_post_test$event_delta.str)
smry$OutputStr  <- paste0(smry$OutputStr,'|vrest_delta|',      ent_v_rest_test$event_delta.str)

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

