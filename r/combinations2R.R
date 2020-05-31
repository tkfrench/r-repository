#added this line
rm(list=ls())
stratVars     <-c('HipKnee', 'Risk', 'Levels', "year") 

stratVarsC   <-rbind(stratVars,  paste0(stratVars,"='All'"))  # -- Add 'All' options to each stratum
#strat.select <-expand.grid(stratVarsC[,1],stratVarsC[,2],stratVarsC[,3],stratVarsC[,4])
  tmp <-sapply(c(1:length(stratVars)),function(x) paste("stratVarsC[,",x,"]",sep=""))
  tmp <-paste0("expand.grid(",paste0(tmp,collapse = ","),")")
strat.select <- eval(parse(text=tmp))    # -- Generate combinations of strata
hold <- subset(strat.select,strat.select[ncol(strat.select)]=='year') # -- Add running 12 mos option
hold[ncol(strat.select)]="year='Running12Mo'"
strat.select <-rbind(strat.select, hold )
strat.select <- apply( strat.select[ ,  ] , 1 , paste , collapse = "," )
strat.select 





#do.call(paste, as.data.frame(strat.select), sep=",")
#paste(strat.select, collapse = ',')
#df <-rbind((as.data.frame(stratVarsC)),(as.data.frame(stratVarsC)))

#df <-rbind(t(as.data.frame(stratVars)),t(as.data.frame(stratVars2)))
#expand.grid(df[,1],df[,2],df[,3],df[,4])

a<-sapply(c(1:length(stratVars)),function(x) paste("df[,",x,"]",sep=""))
#a<-as.character(a)
#a
#b<-paste0("expand.grid(",a,")",sep=",")
#b
#c<-eval(parse(text="expand.grid(df[,1],df[,2],df[,3],df[,4])"))
#eval example
#a <- "5+5"A
#eval(parse(text="5+5"))

#do.call("paste", c(c, sep = ","))
#do.call ("paste",c(t(stratVars),as.list(rep('AAA',16)), sep = ","))
       
#stratVars    <-c('HipKnee',  'Risk',       'Levels', "year")
#lapply(stratVars, function(x) paste0(x, "='All'"))

#stratVarsC <- rbind(stratVars,paste0(stratVars,"=All'"))



