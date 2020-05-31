rm(list=ls())
library(PSJHR)
# following functions will exist in package

tr_OE_v_rest        <-function() {Num    <-.arg1; Exp <-.arg2; alpha <-.arg3
                                  OE_v_rest(Num, Exp, alpha)$all_results_str}

tr_OE_pre_post      <-function() {NumPre <-.arg1; DenPre <-.arg2; ExpPre <-.arg3;
                                  NumPost<-.arg4; DenPost<-.arg5; ExpPost<-.arg6; alpha <-.arg7
                                  OE_pre_post(NumPre, NumPost, ExpPre, ExpPost, alpha)$all_results_str}
tr_welch_t_v_rest   <-function() {Num <-.arg1; Den <-.arg2; NumSqrd <-.arg3; alpha <-.arg4
                                  welch_t_v_rest(Num, Den, NumSqrd,  twosided = TRUE, alpha)$all_results_str}


#testlist =c('OE_test', 'OE_test_PrePost', 'Rate_test', 'Rate_test_PrePost',
#'t_test', 't_test_PrePost', 'NPAR_test')

# some sample data
.arg1 <- c(1,2,3,4,5,6,7,8,9,10)
.arg2 <- c(11,12,13,14,15,16,17,18,19,110); .arg3 <- c(11,12,13,14,15,16,17,18,19,110)
.arg4 <- c(1,2,3,4,5,6,7,8,9,10); .arg5 <- c(11,12,13,14,15,16,17,18,19,110)
.arg6 <- c(1,2,3,4,5,6,7,8,9,10); .arg7 <- c(11,12,13,14,15,16,17,18,19,110)
.arg8 <- rep(0.05, 10)

# Examples:
# Tableau Script call: SCRIPT_STR('tr_OE_v_rest()',   sum([Mort_num]), count([Mort_num]), sum([Mort_exp]), [alpha])
tr_OE_v_rest()
# Tableau Script call: SCRIPT_STR('tr_OE_pre_post()', sum([Mort_num_pre]),  sum([Mort_den_pre], sum([Mort_exp_pre]), 
#                                                     sum([Mort_num_prost]),sum([Mort_den_post],sum([Mort_exp_pre]), [alpha]) 
tr_OE_pre_post()

tr_welch_t_v_rest()