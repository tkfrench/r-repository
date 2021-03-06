logistic simulation1
========================================================
author: Thomas French
date: 12/30/2018
autosize: true

First Slide
========================================================

For more details on authoring R presentations please visit <https://support.rstudio.com/hc/en-us/articles/200486468>.

 Routine for generating sample random dataset for testing Logistic model  
 Generates two variables - y: response variable; x: risk factor  
 q: ~Prevalence of y (set by user)  
 p: ~Prevalence of x (set by user)  
#  rho: Spearman rank correlation
#       used as strength of association of y with x, (set by user).
#       Higher the abs(rho) => higher ROC for glm(y ~ x,family = binomial(link='logit'))
# Four hospital sets of 1,000 cases generated with varying leakage rates% (0,10,20,40)
##################################################################
- Bullet 1
- Bullet 2
- Bullet 3

Slide With Code
========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](logistic simulation1-figure/unnamed-chunk-2-1.png)
