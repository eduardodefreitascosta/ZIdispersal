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


wd <- getwd()

set.seed(13)

source(here("Scripts","function_Weibull_No_p1.r"),local=TRUE)


modelos <- c("dados_No_p1", "modelos_No_p1")
N<-c(200,400,600,800,1000)


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



dir.create(paste(wd,"/Figure_weibull",sep=""))

###-------------------------------------------------------------------------------------------------
# Define files

modelos <- c("dados_No_p1", "modelos_No_p1")
N<-c(200,400,600,800,1000)


###--------------------------------------------------------------------------------------------###
# Datasets generated 
###--------------------------------------------------------------------------------------------###
#Set the covariates.
#Set the first parameter scenario.
# beta1 = 0.5 
# beta2 = 0.5 
# beta3 = 1.5
# beta4 =   2
# beta5 =  -3
# beta6 =   1
#beta7 =  -2
#beta8 = 0.75
#varp=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)
# varp=c(beta1,beta2,beta3,beta4,beta5,beta6)
lime = 10

N<-c(200,400,600,800,1000)


##Remember the parameters

beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1
varp=c(beta1,beta2,beta3,beta4,beta5,beta6)

##Graphics

samp<-c(200,400,600,800,1000)

#Params

betas200 <- read.table(paste(wd,'/','/Output/simulation',N[1],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas400 <- read.table(paste(wd,'/','/Output/simulation',N[2],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas600 <- read.table(paste(wd,'/','/Output/simulation',N[3],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas800 <- read.table(paste(wd,'/','/Output/simulation',N[4],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas1000 <- read.table(paste(wd,'/','/Output/simulation',N[5],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas<-cbind(apply(betas200[,1:6],2,mean),
             apply(betas400[,1:6],2,mean),
             apply(betas600[,1:6],2,mean),
             apply(betas800[,1:6],2,mean),
             apply(betas1000[,1:6],2,mean))

vies<-(betas-varp)

dps<-cbind(apply(betas200[,1:6],2,sd),
           apply(betas400[,1:6],2,sd),
           apply(betas600[,1:6],2,sd),
           apply(betas800[,1:6],2,sd),
           apply(betas1000[,1:6],2,sd))

variancia<-cbind(apply(betas200[,1:6],2,var),
                 apply(betas400[,1:6],2,var),
                 apply(betas600[,1:6],2,var),
                 apply(betas800[,1:6],2,var),
                 apply(betas1000[,1:6],2,var))

var200 <- read.table(paste(wd,"/Output/simulation200/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var400 <- read.table(paste(wd,"/Output/simulation400/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var600 <- read.table(paste(wd,"/Output/simulation600/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var800 <- read.table(paste(wd,"/Output/simulation800/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var1000 <- read.table(paste(wd,"/Output/simulation1000/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))

vars<-cbind(apply(var200[,1:6],2,mean),
            apply(var400[,1:6],2,mean),
            apply(var600[,1:6],2,mean),
            apply(var800[,1:6],2,mean),
            apply(var1000[,1:6],2,mean))

eqm<-(vies^2+variancia)


#Plots

vies1 = t(vies)
vies1 = data.frame(vies1)
attach(vies1)

a <- ggplot(data = vies1, aes(x = samp, y = V1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[1], " )")), x = "Sample Size")

b <- ggplot(data = vies1, aes(x = samp, y = V2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[2], " )")), x = "Sample Size")

c <- ggplot(data = vies1, aes(x = samp, y = V3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[3], " )")), x = "Sample Size")

d <- ggplot(data = vies1, aes(x = samp, y = V4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[4], " )")), x = "Sample Size")

e <- ggplot(data = vies1, aes(x = samp, y = V5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[5], " )")), x = "Sample Size")

f <- ggplot(data = vies1, aes(x = samp, y = V6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[6], " )")), x = "Sample Size")

tiff(file=paste(wd,"/Figure_weibull","/beta_bias.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
gridExtra::grid.arrange(a,b,c,d,e,f,nrow=3)
dev.off()

rownames(vars) = c("Z1","Z2","Z3","Z4","Z5","Z6") 
vars1 = t(vars)
vars1 = data.frame(vars1)
attach(vars1)

a2 <- ggplot(data = vars1, aes(x = samp, y = Z1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[1], " )")), x = "Sample Size")

b2 <- ggplot(data = vars1, aes(x = samp, y = Z2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[2], " )")), x = "Sample Size")

c2 <- ggplot(data = vars1, aes(x = samp, y = Z3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[3], " )")), x = "Sample Size")

d2 <- ggplot(data = vars1, aes(x = samp, y = Z4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[4], " )")), x = "Sample Size")

e2 <- ggplot(data = vars1, aes(x = samp, y = Z5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[5], " )")), x = "Sample Size")

f2 <- ggplot(data = vars1, aes(x = samp, y = Z6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[6], " )")), x = "Sample Size")

tiff(file=paste(wd,"/Figure_weibull","/beta_var.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
gridExtra::grid.arrange(a2,b2,c2,d2,e2,f2,nrow=3)
dev.off()


#Plots

rownames(betas) = c("W1","W2","W3","W4","W5","W6") 
betas1 = t(betas)
betas1 = data.frame(betas1)
attach(betas1)

a3 <- ggplot(data = betas1, aes(x = samp, y = W1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta1) +
  labs(y = expression(paste("Mean (",hat(beta)[1],")" )), x = "Sample Size")

b3 <- ggplot(data = betas1, aes(x = samp, y = W2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta2) +
  labs(y = expression(paste("Mean (",hat(beta)[2],")" )), x = "Sample Size")

c3 <- ggplot(data = betas1, aes(x = samp, y = W3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta3) +
  labs(y = expression(paste("Mean (",hat(beta)[3],")" )), x = "Sample Size")

d3 <- ggplot(data = betas1, aes(x = samp, y = W4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta4) +
  labs(y = expression(paste("Mean (",hat(beta)[4],")" )), x = "Sample Size")

e3 <- ggplot(data = betas1, aes(x = samp, y = W5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta5) +
  labs(y = expression(paste("Mean (",hat(beta)[5],")" )), x = "Sample Size")

f3 <- ggplot(data = betas1, aes(x = samp, y = W6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta6) +
  labs(y = expression(paste("Mean (",hat(beta)[6],")" )), x = "Sample Size")


tiff(file=paste(wd,"/Figure_weibull","/beta_mean.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
grid.arrange(arrangeGrob(a3,b3, ncol=2), 
             arrangeGrob(c3,d3, ncol=2),
             arrangeGrob(e3,f3, ncol=2), ncol=1)

dev.off()


#########################################################################

dir.create(paste(wd,"/Tables_weibull",sep=""))

#Theoretic E(t|x=)

media<-function(x){
  alfa=exp(beta1+beta2*x)
  teta=exp(beta3+beta4*x)
  p0=exp(beta5+beta6*x)/(1+exp(beta5+beta6*x))
  
  media.1=(1-p0)*teta*gamma(1+(1/alfa))
  return(rbind(media.1,p0))
}
media(x=-1)
media(x=1)

#Mean time ZIWeibull

mediaw<-function(x){
  alfa1=numeric()
  teta1=numeric()
  po1=numeric()
  
  for (j in 1:5){
    alfa1[j]=exp(betas[1,j]+betas[2,j]*x) #k in the article
    teta1[j]=exp(betas[3,j]+betas[4,j]*x) #lambda in the article
    po1[j]=exp(betas[5,j]+betas[6,j]*x)/(1+exp(betas[5,j]+betas[6,j]*x))
  }
  
  media1=(1-po1)*teta1*gamma(1+(1/alfa1))
  return (rbind(media1,po1))
}

means<-rbind(
mediaw(x=1),
mediaw(x=-1))

colnames(means)[1:5]<-c("200","400","600","800","1000")

row.names(means)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(means,here("Tables_weibull","means.txt"),sep=";",row.names = T)

###Tables

##Weibull part

##Cover probability

betas200e <- betas200[,1:6]
betas400e <- betas400[,1:6]
betas600e <- betas600[,1:6]
betas800e <- betas800[,1:6]
betas1000e <- betas1000[,1:6]

LimInf200 <- betas200e[,]-1.96*sqrt(variancia[,1])
LimSup200 <- betas200e[,]+1.96*sqrt(variancia[,1])
R200 <- nrow(betas200e)
lim_inf_r200 <- matrix(NA,R200,6)
lim_sup_r200 <- matrix(NA,R200,6)

for (i in 1:6){
  lim_inf_r200[,i] <- betas200e[,i] -1.96*sd(betas200e[,i])
  lim_sup_r200[,i] <- betas200e[,i] +1.96*sd(betas200e[,i])
}
Prob_cobertura200 <- NULL
for (i in 1:6){
  Prob_cobertura200[i] <- (sum(lim_sup_r200[,i] >= varp[i] & lim_inf_r200[,i] <= varp[i]))/R200
}


LimInf400 <- betas400e[,]-1.96*sqrt(variancia[,2])
LimSup400 <- betas400e[,]+1.96*sqrt(variancia[,2])
R400 <- nrow(betas400e)
lim_inf_r400 <- matrix(NA,R400,6)
lim_sup_r400 <- matrix(NA,R400,6)

for (i in 1:6){
  lim_inf_r400[,i] <- betas400e[,i] -1.96*sd(betas400e[,i])
  lim_sup_r400[,i] <- betas400e[,i] +1.96*sd(betas400e[,i])
}
Prob_cobertura400 <- NULL
for (i in 1:6){
  Prob_cobertura400[i] <- (sum(lim_sup_r400[,i] >= varp[i] & lim_inf_r400[,i] <= varp[i]))/R400
}


LimInf600 <- betas600e[,]-1.96*sqrt(variancia[,3])
LimSup600 <- betas600e[,]+1.96*sqrt(variancia[,3])
R600 <- nrow(betas600e)
lim_inf_r600 <- matrix(NA,R600,6)
lim_sup_r600 <- matrix(NA,R600,6)

for (i in 1:6){
  lim_inf_r600[,i] <- betas600e[,i] -1.96*sd(betas600e[,i])
  lim_sup_r600[,i] <- betas600e[,i] +1.96*sd(betas600e[,i])
}
Prob_cobertura600 <- NULL
for (i in 1:6){
  Prob_cobertura600[i] <- (sum(lim_sup_r600[,i] >= varp[i] & lim_inf_r600[,i] <= varp[i]))/R600
}


LimInf800 <- betas800e[,]-1.96*sqrt(variancia[,4])
LimSup800 <- betas800e[,]+1.96*sqrt(variancia[,4])
R800 <- nrow(betas800e)
lim_inf_r800 <- matrix(NA,R800,6)
lim_sup_r800 <- matrix(NA,R800,6)

for (i in 1:6){
  lim_inf_r800[,i] <- betas800e[,i] -1.96*sd(betas800e[,i])
  lim_sup_r800[,i] <- betas800e[,i] +1.96*sd(betas800e[,i])
}
Prob_cobertura800 <- NULL
for (i in 1:6){
  Prob_cobertura800[i] <- (sum(lim_sup_r800[,i] >= varp[i] & lim_inf_r800[,i] <= varp[i]))/R800
}

LimInf1000 <- betas1000e[,]-1.96*sqrt(variancia[,5])
LimSup1000 <- betas1000e[,]+1.96*sqrt(variancia[,5])
R1000 <- nrow(betas1000e)
lim_inf_r1000 <- matrix(NA,R1000,6)
lim_sup_r1000 <- matrix(NA,R1000,6)

for (i in 1:6){
  lim_inf_r1000[,i] <- betas1000e[,i] -1.96*sd(betas1000e[,i])
  lim_sup_r1000[,i] <- betas1000e[,i] +1.96*sd(betas1000e[,i])
}
Prob_cobertura1000 <- NULL
for (i in 1:6){
  Prob_cobertura1000[i] <- (sum(lim_sup_r1000[,i] >= varp[i] & lim_inf_r1000[,i] <= varp[i]))/R1000
}

probcob <- cbind(Prob_cobertura200,
                 Prob_cobertura400,
                 Prob_cobertura600,
                 Prob_cobertura800,
                 Prob_cobertura1000)


tab200<-cbind(betas[,1],cbind(betas[,1]-1.96*sqrt(variancia[,1]),betas[,1]+1.96*sqrt(variancia[,1])),dps[,1],sqrt(vars[,1]),vies[,1],eqm[,1],probcob[,1])
row.names(tab200)<-c(1,2,3,4,5,6)

tab400<-cbind(betas[,2],cbind(betas[,2]-1.96*sqrt(variancia[,2]),betas[,2]+1.96*sqrt(variancia[,2])),dps[,2],sqrt(vars[,2]),vies[,2],eqm[,2],probcob[,2])
row.names(tab400)<-c(1,2,3,4,5,6)

tab600<-cbind(betas[,3],cbind(betas[,3]-1.96*sqrt(variancia[,3]),betas[,3]+1.96*sqrt(variancia[,3])),dps[,3],sqrt(vars[,3]),vies[,3],eqm[,3],probcob[,3])
row.names(tab600)<-c(1,2,3,4,5,6)

tab800<-cbind(betas[,4],cbind(betas[,4]-1.96*sqrt(variancia[,4]),betas[,4]+1.96*sqrt(variancia[,4])),dps[,4],sqrt(vars[,4]),vies[,4],eqm[,4],probcob[,4])
row.names(tab800)<-c(1,2,3,4,5,6)

tab1000<-cbind(betas[,5],cbind(betas[,5]-1.96*sqrt(variancia[,5]),betas[,5]+1.96*sqrt(variancia[,5])),dps[,5],sqrt(vars[,5]),vies[,5],eqm[,5],probcob[,5])
row.names(tab1000)<-c(1,2,3,4,5,6)



final1<-rbind(tab200,tab400,tab600,tab800,tab1000)
row.names(final1)<-rep(c(1:6),5)
colnames(final1)<-c("media","LCI","UCI","SD of betas","Mean Std Err","vies","eqmg","CP")

write.table(final1,here("Tables_weibull","final.txt"),sep=";",row.names = F)

