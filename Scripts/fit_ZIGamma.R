rm(list=ls(all=TRUE)) 

#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","flexsurv","knitr","glmmsr","plotly","gridExtra","grid","ggridges","ggthemes","summarytools","ggcorrplot","survival")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


wd <- getwd()

set.seed(13)
beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1


###--------------------------------------------------------------------------------------------###
# #Gamma model (weibull data) 
###--------------------------------------------------------------------------------------------###

nomes<-c("linear","gamma","weibull")
replica<-1000
N<-c(200,400,600,800,1000)

for (j in 1:length(N)){
  cat("\n\n Sample: ", N[j], "\n")
  
  wd.sim<-paste(wd,"/Output/simulation",N[j],"/modelos_No_p1/", 'gamma', sep="")  
  param_Gamma<- file.path(wd.sim, "param_Gamma.txt") 
  Erro_Gamma <- file.path(wd.sim, "Erro_Gamma.txt")
  
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
    
    custom.gamma <- list(name = "gamma",
                         pars = c("shape", "scale"),
                         location = "scale",
                         transforms = c(log, log),
                         inv.transforms = c(exp, exp),
                         inits = function(t){
                           c(0.01, 100)},
                         method="Nelder-Mead"
    )
    
    c<-flexsurvreg(Surv(t, d,type='right')~x,
                   anc = list(shape = ~ x),data=dados,dist=custom.gamma,
                   subset=dados$z==0,
                   control=list(fnscale = 2500)
    )
    
    
    varic=try(solve(c$opt$hessian))
    
    param_gamma<- file.path(wd.sim, "param_gamma.txt") 
    
    Erro_gamma <- file.path(wd.sim, "Erro_gamma.txt")
    
    cat(c(c$res.t[1],c$res.t[4],c$res.t[2],c$res.t[3],a$coefficients,c$loglik,c$AIC),file = param_gamma,append = TRUE, "\n") 
    cat(c(varic[1,1],varic[4,4],varic[2,2],varic[3,3],b$coefficients[1,2]^2,b$coefficients[2,2]^2),file = Erro_gamma,append = TRUE, "\n")   }
  
}

###--------------------------------------------------------------------------------------------###
# #Gamma model (gamma data) 
###--------------------------------------------------------------------------------------------###

replica<-1000
N<-c(200,400,600,800,1000)
nomes<-c("linear","gamma","weibull")
for (j in 1:length(N)){
  cat("\n\n Sample: ", N[j], "\n")
  
  dir.create( paste(wd, "/Output/simulation_gamma",N[j],"/modelos_No_p1/", sep="") )
  dir.create( paste(wd, "/Output/simulation_gamma",N[j],"/modelos_No_p1/","gamma", sep="") )
  wd.sim<-paste(wd,"/Output/simulation_gamma",N[j],"/modelos_No_p1/", "gamma", sep="")
  
  
  for(iteracao in 1:replica){ 
    
    tempo_inicial = proc.time()
    cat("\n\n iteraction: ", iteracao, "\n")  
    data.labels <- paste("dados",1:replica, ".txt", sep="")   
    
    wd.dados <- paste(wd, "/Output/simulation_gamma",N[j],"/dados_No_p1",sep="")
    data.name <- data.labels[iteracao]
    data.local  <- file.path(wd.dados,data.name)
    dados <- read.table(data.local, head=TRUE)
    
    
    a<-glm(dados$z ~ dados$x,family=binomial(link='logit'))
    b<-summary(a)
    
    custom.gamma <- list(name = "gamma",
                           pars = c("shape", "scale"),
                           location = "scale",
                           transforms = c(log, log),
                           inv.transforms = c(exp, exp),
                           inits = function(t){
                             c(0.01, 100)},
                         method="Nelder-Mead"
                          )
    
    c<-flexsurvreg(Surv(t, d,type='right')~x,
                   anc = list(shape = ~ x),data=dados,dist=custom.gamma,
                   subset=dados$z==0,
                   control=list(fnscale = 2500)
                   )
    
    
    varic=try(solve(c$opt$hessian))
      
      param_gamma<- file.path(wd.sim, "param_gamma.txt") 
      
      Erro_gamma <- file.path(wd.sim, "Erro_gamma.txt")
  
    cat(c(c$res.t[1],c$res.t[4],c$res.t[2],c$res.t[3],a$coefficients,c$loglik,c$AIC),file = param_gamma,append = TRUE, "\n") 
    cat(c(varic[1,1],varic[4,4],varic[2,2],varic[3,3],b$coefficients[1,2]^2,b$coefficients[2,2]^2),file = Erro_gamma,append = TRUE, "\n") 
  }
  
}


