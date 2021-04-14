rm(list=ls(all=TRUE)) 

#Packages to be used
packages<-c("readxl","here","knitr","tidyverse","ggplot2","flexsurv","knitr","glmmsr","plotly","gridExtra","grid","ggridges","ggthemes","summarytools","ggcorrplot")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Set location and seed
wd <- getwd()

set.seed(13)


###--------------------------------------------------------------------------------------------###
# Set the function for generation 
###--------------------------------------------------------------------------------------------###

model <- function(N,varp){
  a       =  exp(varp[1])*exp(y*varp[2])
  lam     =  exp(varp[3])*exp(y*varp[4])
  
  # Parameters - WEIGHT IN ZEROS AND INFS: GAMMA ZERO AND GAMMA ONES WITH A COVARIABLE
  # gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6]))  
  # gamm_inf = (exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  
  U=runif(N)
  V=runif(N,gamm_zer,1)
  y=numeric(N)
  for(i in 1:N){
    if(U[i]<=gamm_zer[i]){ 
      y[i] = 0  
    }else 
      if(U[i]>gamm_zer[i] & U[i] <= 1){
        #y[i]=qweibull((V[i]-gamm_zer[i])/(1-gamm_zer[i]), a[i], lam[i])
        y[i] = lam[i]*((-log(1 - (V[i]-gamm_zer[i])/(1-gamm_zer[i])))^(1/a[i]))
      }else {      
        y[i]=+Inf
      }
  }
  return(y)
}



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


varp=c(beta1,beta2,beta3,beta4,beta5,beta6)
lime = 10
N<-c(200,400,600,800,1000)
replica <- 1000

tempo_inicial = proc.time()

#set.seed(20)
for (j in 1:length(N)){
  tempo_inicial = proc.time()
  
  for(iteracao in 1:replica){
    y=rnorm(N[j],0,1)
    tp = model(N[j],varp)
    tM = max(tp[tp<+Inf])
  
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

