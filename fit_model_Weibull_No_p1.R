###--------------------------------------------------------------------------------------------###
# The estimation steps follow the sequence outlined in the paper about dependent censoring, 
# for Weibull adjustment.
###--------------------------------------------------------------------------------------------###
rm(list=ls(all=TRUE)) 


if (!require(flexsurv)){
  install.packages("flexsurv")
}
library(flexsurv)


wd <- getwd()

set.seed(13)

source("function_Weibull_No_p1.r",local=TRUE)

###-------------------------------------------------------------------------------------------------
# Define files

modelos <- c("dados_No_p1", "modelos_No_p1")
N<-c(200,400,600,800,1000)

dir.create(paste(wd,"/Output",sep=""))

  
  for(j in 1:length(N)){
    dir.create( paste(wd, "/Output/simulation",N[j],"/", sep="") )
    
  } 

for(j in 1:length(N)){
for (i in 1:length(modelos)){
  
  dir.create( paste(wd, "/Output/simulation",N[j],"/", modelos[i], sep="") )  
}
}



###--------------------------------------------------------------------------------------------###
# Datasets generated 
###--------------------------------------------------------------------------------------------###
#Set the covariates.
#Set the first parameter scenario.
beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1
#beta7 =  -2
#beta8 = 0.75
#varp=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)
varp=c(beta1,beta2,beta3,beta4,beta5,beta6)
lime = 10
N<-c(200,400,600,800,1000)
replica <- 1000

tempo_inicial = proc.time()

set.seed(13)
for (j in 1:length(N)){
  tempo_inicial = proc.time()
  
  for(iteracao in 1:replica){
  y=rnorm(N[j],0,1)
  tp = model(N[j],varp)
  tM = max(tp[tp<+Inf])
  #Z=runif(N[j],0,tM) #time until failure
  Z=rgamma(N[j], shape=0.25, scale = 128.2) #time until censor
  
  time=numeric(); delta=numeric(); gamma=numeric()
  
  for (i in 1:N[j]){
    time[i] = min(tp[i],Z[i]) #observed time (minimum between failure time and censor)
  
    
    if (time[i] < Z[i]){
      delta[i]=1          #delta: indicative for failure and censor
    } else {delta[i]=0}
    
    if (time[i]==0){
      gamma[i]=1
    }else{gamma[i]=0}
  
    }
    
  t<-time  #Observed time;
  d<-delta #Indicative for failure and cens
  x<-y     #covariable
  z<-gamma #Zero indicative
  
  dados <- cbind(t,d,z,x)
  
  data.labels <- paste("dados",1:replica, ".txt", sep="")
 
  wd.dados <- paste(wd, "/Output/simulation",N[j],"/dados_No_p1",sep="")
  
  data.local  <- file.path(wd.dados, data.labels[iteracao])
 
  write.table(dados, file = data.local, row.names=FALSE, col.names=TRUE, append=FALSE)
 }
}




###--------------------------------------------------------------------------------------------###
# Monte Carlo Study
###--------------------------------------------------------------------------------------------###
nomes<-c("linear","gamma","weibull")
N<-c(200,400,600,800,1000)
replica <- 1000
continue <- TRUE

for (j in 1:5){
cat("\n\n Sample: ", N[j], "\n")
 for(iteracao in 1:replica){ 
	
   tempo_inicial = proc.time()
   cat("\n\n Iteraction: ", iteracao, "\n")  
   data.labels <- paste("dados",1:replica, ".txt", sep="")   
   
   wd.dados <- paste(wd, "/Output/simulation",N[j],"/dados_No_p1",sep="")
   data.name <- data.labels[iteracao]
   data.local  <- file.path(wd.dados,data.name)
   dados <- read.table(data.local, head=TRUE)
   
   time<-dados$t
   delta<-dados$d
   y<-dados$x
   
   #varp=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)
   varp=c(beta1,beta2,beta3,beta4,beta5,beta6)
   
   fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))
   estc = try(fit$par)
   Hc=try(fit$hessian); varic=try(-solve(Hc))
 
 
     for(i in 1:length(nomes)){
     
       dir.create( paste(wd, "/Output/simulation",N[j],"/modelos_No_p1/", nomes[i], sep="") )
       wd.sim<-paste(wd,"/Output/simulation",N[j],"/modelos_No_p1/", nomes[i], sep="")
       
       
       param_Weibull<- file.path(wd.sim[1], "param_Weibull.txt") 
       
       Erro_Weibull <- file.path(wd.sim[1], "Erro_Weibull.txt")
       
     }
 
 
 
    if ( is.finite(estc[1]) & is.finite(estc[2]) & is.finite(estc[3]) & is.finite(estc[4]) & 
          is.finite(estc[5]) & is.finite(estc[6])) {
        if ( is.finite(varic[1,1]) & is.finite(varic[2,2]) & is.finite(varic[3,3]) & is.finite(varic[4,4]) & 
            is.finite(varic[5,5]) & is.finite(varic[6,6]) ){
          if ( varic[1,1] > 0 & varic[2,2] > 0 & varic[3,3] > 0 & varic[4,4] > 0 & 
                varic[5,5] > 0 & varic[6,6] > 0  &
                estc[1] < lime & estc[2] < lime & estc[3] < lime & estc[4] < lime &
                estc[5] < lime & estc[6] < lime ){
       
       var.varp1 =varic[1,1]; Lvarp1 =estc[1]-1.96*sqrt(var.varp1); Uvarp1=estc[1]+1.96*sqrt(var.varp1)
       var.varp2 =varic[2,2]; Lvarp2 =estc[2]-1.96*sqrt(var.varp2); Uvarp2=estc[2]+1.96*sqrt(var.varp2)
       var.varp3 =varic[3,3]; Lvarp3 =estc[3]-1.96*sqrt(var.varp3); Uvarp3=estc[3]+1.96*sqrt(var.varp3)
       var.varp4 =varic[4,4]; Lvarp4 =estc[4]-1.96*sqrt(var.varp4); Uvarp4=estc[4]+1.96*sqrt(var.varp4)
       var.varp5 =varic[5,5]; Lvarp5 =estc[5]-1.96*sqrt(var.varp5); Uvarp5=estc[5]+1.96*sqrt(var.varp5)
       var.varp6 =varic[6,6]; Lvarp6 =estc[6]-1.96*sqrt(var.varp6); Uvarp6=estc[6]+1.96*sqrt(var.varp6)
       #var.varp7 =varic[7,7]; Lvarp7 =estc[7]-1.96*sqrt(var.varp7); Uvarp7=estc[7]+1.96*sqrt(var.varp7)
       #var.varp8 =varic[8,8]; Lvarp8 =estc[8]-1.96*sqrt(var.varp8); Uvarp8=estc[8]+1.96*sqrt(var.varp8)
       
       cat(c(estc,fit$value),file = param_Weibull,append = TRUE, "\n") 
       cat(c(var.varp1,var.varp2,var.varp3,var.varp4,var.varp5,var.varp6),file = Erro_Weibull,append = TRUE, "\n") 
       
       }
     }  
    }
 
  } 
}



##Gamma model

replica<-1000
N<-c(200,400,600,800,1000)

for (j in 1:length(N)){
  cat("\n\n Sample: ", N[j], "\n")
  for(iteracao in 1:replica){ 
    
    tempo_inicial = proc.time()
    cat("\n\n iteraction: ", iteracao, "\n")  
    data.labels <- paste("dados",1:replica, ".txt", sep="")   
    
    wd.dados <- paste(wd, "/Output/simulation",N[j],"/dados_No_p1",sep="")
    data.name <- data.labels[iteracao]
    data.local  <- file.path(wd.dados,data.name)
    dados <- read.table(data.local, head=TRUE)
    
    
    
      a<-glm(dados$z ~ dados$x,family=binomial(link='logit'))
      b<-summary(a)
      c<-flexsurvreg(Surv(dados$t, dados$d,type='right')~dados$x,dist="gamma",subset=dados$z==0, method="Nelder-Mead",control=list(fnscale = 250000),integ.opts = list(rel.tol=1e-8))
      
      varic=try(solve(c$opt$hessian))
      
      wd.sim<-paste(wd,"/Output/simulation",N[j],"/modelos_No_p1/", 'gamma', sep="")  
      param_Gamma<- file.path(wd.sim, "param_Gamma.txt") 
      Erro_Gamma <- file.path(wd.sim, "Erro_Gamma.txt")
      
      cat(c(a$coefficients,c$res[1],c$res[2],c$res[3], c$loglik,c$AIC),file = param_Gamma,append = TRUE, "\n") 
      cat(c(b$coefficients[1,2],b$coefficients[2,2],varic[1,1],varic[2,2],varic[3,3]),file = Erro_Gamma,append = TRUE, "\n") 
  }

}

source("fit_model_Weibull_No_p1_graficos_ggplot2.r",local=TRUE)
