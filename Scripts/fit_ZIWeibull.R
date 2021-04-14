###--------------------------------------------------------------------------------------------###
# The estimation steps follow the sequence outlined in the paper about dependent censoring, 
# for Weibull adjustment.
###--------------------------------------------------------------------------------------------###

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



###--------------------------------------------------------------------------------------------###
# Monte Carlo Study (Weibull data)
###--------------------------------------------------------------------------------------------###
nomes<-c("linear","gamma","weibull")
N<-c(200,400,600,800,1000)
modelos <- c("dados_No_p1", "modelos_No_p1")
replica <- 1000
continue <- TRUE

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
tempo_inicial = proc.time()

wd <- getwd()

set.seed(13)

source(here("Scripts","function_Weibull_No_p1.r"),local=TRUE)

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




###--------------------------------------------------------------------------------------------###
# Monte Carlo Study (gamma data)
###--------------------------------------------------------------------------------------------###
nomes<-c("linear","gamma","weibull")
N<-c(200,400,600,800,1000)
modelos <- c("dados_No_p1", "modelos_No_p1")
replica <- 1000
continue <- TRUE

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
tempo_inicial = proc.time()

wd <- getwd()

set.seed(13)

source(here("Scripts","function_Weibull_No_p1.r"),local=TRUE)


for (j in 1:5){
  cat("\n\n Sample: ", N[j], "\n")
  for(iteracao in 1:replica){ 
    
    tempo_inicial = proc.time()
    cat("\n\n Iteraction: ", iteracao, "\n")  
    data.labels <- paste("dados",1:replica, ".txt", sep="")   
    
    wd.dados <- paste(wd, "/Output/simulation_gamma",N[j],"/dados_No_p1",sep="")
    data.name <- data.labels[iteracao]
    data.local  <- file.path(wd.dados,data.name)
    dados <- read.table(data.local, head=TRUE)
    
    time<-dados$t
    delta<-dados$d
    y<-dados$x
    
    varp=c(beta1,beta2,beta3,beta4,beta5,beta6)
    
    fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))
    estc = try(fit$par)
    Hc=try(fit$hessian); varic=try(-solve(Hc))
    
    
    for(i in 1:length(nomes)){
      
      dir.create( paste(wd, "/Output/simulation_gamma",N[j],"/modelos_No_p1/", nomes[i], sep="") )
      wd.sim<-paste(wd,"/Output/simulation_gamma",N[j],"/modelos_No_p1/", nomes[i], sep="")
      
      
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




