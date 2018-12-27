rm(list=ls())
library(plyr)

stratumVars  <- list("gender", "age_gt65", 'year', 'fracture')
dontCareVars <- sapply(stratumVars, paste, "DK", sep='')

listIn<-list()
for(i in 1:length(stratumVars)) {
  listIn[[length(listIn)+1]] <- combn(stratumVars,i)
  #print(listIn)
}

varCombos<-list()
for (j in 1:length(listIn)){
  for (k in 1:length(listIn[[j]][1,])){
    # print(list(listIn[[j]][,k]))
   #varCombos[[length(varCombos)+1]] <- list(listIn[[j]][,k])  
    varCombos[[length(varCombos)+1]] <- list(c(unlist(listIn[[j]][,k])))
    #varCombos[[length(varCombos)]] <- list(c(unlist(varCombos[[length(varCombos)]]),unlist(dontCareVars)))
  }
}
varCombos
for (j in 1:length(varCombos)){
  for (k in 1:length(stratumVars){
    
    if(!grepl(stratumVars[k],varCombos[j])){
      varCombos[j]<-append((varCombos[j]),paste0(stratumVars[k],'ALL'))
    }
  }
}


a<-varCombos[15]
 stratumVars

listIn
length(listIn)        # how many lists
length(listIn[[1]][1,]) # how many varComboss in list 1
list(listIn[[1]][,1]) #List 1,1
list(listIn[[1]][,2]) #List 1,2
list(listIn[[1]][,3]) #List 1,3
list(listIn[[1]][,4]) #List 1,4
length(listIn[[2]][1,]) # how many varComboss in list 2
list(listIn[[2]][,1]) #List 2,1
list(listIn[[2]][,2]) #List 2,2
list(listIn[[2]][,3]) #List 2,3
list(listIn[[2]][,4]) #List 2,4
list(listIn[[2]][,5]) #List 2,5
list(listIn[[2]][,6]) #List 2,6
length(listIn[[3]][1,]) # how many varComboss in list 3
list(listIn[[3]][,1]) #List 3,1
list(listIn[[3]][,2]) #List 3,2
list(listIn[[3]][,3]) #List 3,3
list(listIn[[3]][,4]) #List 3,4
length(listIn[[4]][1,]) # how many varComboss in list 4
list(listIn[[4]][,1]) #List 4,1

varCombos<-list()
for (j in 1:length(listIn)){
  for (k in 1:length(listIn[[j]][1,])){
   # print(list(listIn[[j]][,k]))
    varCombos[[length(varCombos)+1]] <- list(listIn[[j]][,k])   
  }
}