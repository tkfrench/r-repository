rm(list=ls())
require(RODBC)
require(exact2x2) # Revised Fisher exact test
require(plyr) 
require(dplyr) 
#====================================
# Generate sample cohort
# SQL Query goes here
# ===================================
rndmnum =sample(1:100,10000,replace=T)
mydata <- data.frame(
  # Entities--------------
  entityType =paste0("",length(rndmnum)),
  entityID   =paste0("",length(rndmnum)),
  provider=paste0("p",1+ rndmnum %%100),     # calssify among 100 providers    
  ministry=paste0("h",1+ rndmnum %%24),      # calssify among  24 Ministries
  region  =paste0("R",1+ rndmnum %%5),       # calssify among   5 Regions
  # year-----------------
  year    =sample(c(2016,2017,2018),length(rndmnum),replace=T), 
  # Stratum--------------
  hip_knee=sample(c('Hip','Knee'),length(rndmnum),replace=T), 
  fracture=sample(c(0,1),length(rndmnum),replace=T) ,
  # Outcomes-------------- 
  binaryoutome  =sample(c(0,1),length(rndmnum),replace=T), 
  intervaloutome=sample(c(0,1),length(rndmnum),replace=T) 
)
rm(rndmnum)
# add a missing value
mydata$binaryoutome[1]=NA
head(mydata)

#===============================
#   Define a few functions here
#===============================
# - fisher exact test --
# --  set muex= FALSE to compare vs System
Bprob <- function(a,b,c,d,muex=TRUE){
    if(!muex){c=c-a; d=d-b}
    round(binomMeld.test(a,b,c,d)$p.value, 5)}
Bprob(5,10,25,100,muex=T)
Bprob(5,10,25,100,muex=F)

# -Welch test --
welchsmryttest <- function(m1,m2,sd1,sd2,num1,num2,twosided=TRUE,alpha=0.05){
  se <- sqrt(sd1*sd1/num1+sd2*sd2/num2)
  t  <- (m1-m2)/se
  df <-(sd1*sd1/num1+sd2*sd2/num2)^2 / (sd1^4/(num1^2*(num1-1))  + sd2^4/(num2^2*(num2-1)))
  if (twosided) {
    pvalue=2*(1-pt(abs(t),df=df))
  } else {
    pvalue=   1-pt(abs(t),df=df)  
  }  
  status <- ifelse(pvalue <= alpha ,ifelse((m1-m2) <0, "Lower","Higher"),"NS")
  #return <- list(pvalue=pvalue, t=t, df=df, status=status,alpha)  
  return(round(pvalue,digits=5))
}

#=============================================
# Build/append records to capture all combinations
# including collapsed startum (e.g, All months and All fractures, etc.)
#=============================================
df<- mydata
t<-mydata; t$year<-"All";                                        df <-rbind(t, df)
t<-mydata;                 t$hip_knee<-"All";                    df <-rbind(t, df)
t<-mydata;                                    t$fracture<-"All"; df <-rbind(t, df)
t<-mydata; t$year<-"All";  t$hip_knee<-"All";                    df <-rbind(t, df)
t<-mydata; t$year<-"All";                     t$fracture<-"All"; df <-rbind(t, df)
t<-mydata;                 t$hip_knee<-"All"; t$fracture<-"All"; df <-rbind(t, df)
t<-mydata; t$year<-"All";  t$hip_knee<-"All"; t$fracture<-"All"; df <-rbind(t, df)
rm(t)
#=============================================
# Generate aggregate values for each outcome
#   by entity type and entity ids
#=============================================


outcomegroup <- c("binaryoutome","intervaloutome")
stratumgroup <- c("year","hip_knee", "fracture")
# System --
strata <- c(stratumgroup)
sys_stats <- df %>%
  group_by_at(.vars=strata)%>%
  summarise_at(vars(outcomegroup), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE)
sys_stats$entityType <-'System'
sys_stats$entityID   <-'System'
sys_stats$region     <-'NA'
sys_stats$ministry   <-'NA'
sys_stats$provider   <-'NA'
sys_stats <- sys_stats %>% select(entityType, region, ministry, provider, everything()) #reorder
head(sys_stats)

# Region --
strata <- c(stratumgroup,"region")
rgn_stats <- df %>%
  group_by_at(.vars=strata)%>%
  summarise_at(vars(outcomegroup), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE)
rgn_stats$entityType <-'Region'
rgn_stats$entityID   <- rgn_stats$region
rgn_stats$ministry   <-'NA'
rgn_stats$provider   <-'NA'
rgn_stats <- rgn_stats %>% select(entityType, region, ministry, provider, everything()) #reorder
head(rgn_stats)

# Ministry --
strata <- c(stratumgroup,"region","ministry")
min_stats <- df %>%
  group_by_at(.vars=strata)%>%
  summarise_at(vars(outcomegroup), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE)
min_stats$entityType <-'ministry'
min_stats$entityID   <- min_stats$region
min_stats$ministry   <- min_stats$ministry
min_stats$provider   <-'NA'
min_stats <- min_stats %>% select(entityType, region, ministry, provider, everything()) #reorder
head(min_stats)

# Provider --
strata <- c(stratumgroup,"region","ministry","provider")
prv_stats <- df %>%
  group_by_at(.vars=strata)%>%
  summarise_at(vars(outcomegroup), funs(mean,sum,sd, n=sum(!is.na(.))), na.rm = TRUE)
prv_stats$entityType <-'provider'
prv_stats$entityID   <- prv_stats$region
prv_stats$ministry   <- prv_stats$ministry
prv_stats$provider   <- prv_stats$provider
prv_stats <- prv_stats %>% select(entityType, region, ministry, provider, everything()) #reorder
head(prv_stats)

#=============================================
# Stack all entity stats 
#  -Also join sys stats to all records
#   to facilitate entity vs system significane tests
#=============================================
s<- subset(sys_stats,select= c(-entityType, -entityID, -region, -ministry, -provider))
names(s) <- c(names(c(s[1:3])), paste0(names(s[4:ncol(s)]),'_sys'))
stats <- rbind(sys_stats, rgn_stats, min_stats, prv_stats)
stats <- merge(stats, s, by = c("year", "hip_knee", "fracture"))
# clean up
rm(mydata, df, s, sys_stats, rgn_stats, min_stats, prv_stats)
head(stats)

#====================================================
# Perform statistical tests on outcome vars vs system
#==================================================== 
# Rates:   Entity vs System
tic=proc.time()[3]
stats$binaryoutome_pvalue <-with(stats, mapply(Bprob, binaryoutome_sum, binaryoutome_n, binaryoutome_sum_sys, binaryoutome_n_sys, muex=F))
toc<-proc.time()[3] - tic; toc

# Averages: Entity vs System
tic=proc.time()[3]
stats$intervaloutome_pvalue <-with(stats, 
 mapply(welchsmryttest, 
 intervaloutome_mean, intervaloutome_mean_sys, 
 intervaloutome_sd,   intervaloutome_sd_sys,
  intervaloutome_n,    intervaloutome_n_sys
  ))
toc<-proc.time()[3] - tic; toc
head(stats)

stratumgroup <- c("year","hip_knee", "fracture")
mystr<-data.frame(ndx=numeric(),x=character(),stringsAsFactors=FALSE)

stratUnits<-lapply(stratumgroup, function(x) paste(x,"=\'All\';",sep=' '))
for(i in 1:2){
 #for(i in 1:length(stratUnits)){  
  a1    <- data.frame(ndx=t(combn(length(stratUnits),i,simplify=T) ))
  a1$x  <- unlist(stratUnits[a1$ndx])
  mystr <- rbind(mystr,a1)
}
a2<-t(combn(length(stratumgroup),2,simplify=T) )
a3<-t(combn(length(stratumgroup),3,simplify=T))
mystr<-data.frame(x=replicate(4,'a'))