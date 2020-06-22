###--------------------------------------------------------------------------------------------###
# The estimation steps follow the sequence outlined in the paper about dependent censoring, 
# for Weibull adjustment.
###--------------------------------------------------------------------------------------------###
rm(list=ls(all=TRUE)) 

#wd <- setwd("C:\\Users\\Usuario\\Dropbox\\Simulacao_nova")
wd <- getwd()

set.seed(123456)


source("function_Weibull_No_p1.r",local=TRUE)

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

betas200 <- read.table(paste(wd,'/','simulation',N[1],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas400 <- read.table(paste(wd,'/','simulation',N[2],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas600 <- read.table(paste(wd,'/','simulation',N[3],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas800 <- read.table(paste(wd,'/','simulation',N[4],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas1000 <- read.table(paste(wd,'/','simulation',N[5],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
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

var200 <- read.table(paste(wd,"/simulation200/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var400 <- read.table(paste(wd,"/simulation400/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var600 <- read.table(paste(wd,"/simulation600/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var800 <- read.table(paste(wd,"/simulation800/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var1000 <- read.table(paste(wd,"/simulation1000/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))

vars<-cbind(apply(var200[,1:6],2,mean),
            apply(var400[,1:6],2,mean),
            apply(var600[,1:6],2,mean),
            apply(var800[,1:6],2,mean),
            apply(var1000[,1:6],2,mean))

eqm<-(vies^2+variancia)


#Plots

library(ggplot2)
library(tidyr)

vies1 = t(vies)
vies1 = data.frame(vies1)
attach(vies1)

a <- ggplot(data = vies1, aes(x = samp, y = V1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[1], " )")), x = "Sample")

b <- ggplot(data = vies1, aes(x = samp, y = V2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[2], " )")), x = "Sample")

c <- ggplot(data = vies1, aes(x = samp, y = V3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[3], " )")), x = "Sample")

d <- ggplot(data = vies1, aes(x = samp, y = V4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[4], " )")), x = "Sample")

e <- ggplot(data = vies1, aes(x = samp, y = V5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[5], " )")), x = "Sample")

f <- ggplot(data = vies1, aes(x = samp, y = V6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[6], " )")), x = "Sample")

gridExtra::grid.arrange(a,b,c,d,e,f,nrow=3)




rownames(vars) = c("Z1","Z2","Z3","Z4","Z5","Z6") 
vars1 = t(vars)
vars1 = data.frame(vars1)
attach(vars1)

a2 <- ggplot(data = vars1, aes(x = samp, y = Z1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[1], " )")), x = "Sample")

b2 <- ggplot(data = vars1, aes(x = samp, y = Z2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[2], " )")), x = "Sample")

c2 <- ggplot(data = vars1, aes(x = samp, y = Z3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[3], " )")), x = "Sample")

d2 <- ggplot(data = vars1, aes(x = samp, y = Z4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[4], " )")), x = "Sample")

e2 <- ggplot(data = vars1, aes(x = samp, y = Z5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[5], " )")), x = "Sample")

f2 <- ggplot(data = vars1, aes(x = samp, y = Z6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[6], " )")), x = "Sample")

gridExtra::grid.arrange(a2,b2,c2,d2,e2,f2,nrow=3)



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


library(gridExtra)
grid.arrange(arrangeGrob(a3,b3, ncol=2), 
                        arrangeGrob(c3,d3, ncol=2),
                        arrangeGrob(e3,f3, ncol=2), ncol=1)




#Modelo gamma

##Graphics
gammas200 <- read.table(paste(wd,'/','simulation',N[1],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas400 <- read.table(paste(wd,'/','simulation',N[2],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas600 <- read.table(paste(wd,'/','simulation',N[3],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas800 <- read.table(paste(wd,'/','simulation',N[4],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas1000 <- read.table(paste(wd,'/','simulation',N[5],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))

#Fazendo as medias dos parametros exponencializados
gammas<-cbind(apply(cbind((gammas200[,3:4]),gammas200[,5]),2,mean),
              apply(cbind((gammas400[,3:4]),gammas400[,5]),2,mean),
              apply(cbind((gammas600[,3:4]),gammas600[,5]),2,mean),
              apply(cbind((gammas800[,3:4]),gammas800[,5]),2,mean),
              apply(cbind((gammas1000[,3:4]),gammas1000[,5]),2,mean))

#Fazendo as  medias dos parametros da parte logistica
logis<-cbind(apply((gammas200[,1:2]),2,mean),
             apply((gammas400[,1:2]),2,mean),
             apply((gammas600[,1:2]),2,mean),
             apply((gammas800[,1:2]),2,mean),
             apply((gammas1000[,1:2]),2,mean))

#Fazendo os SD dos parametros exponencializados
gdps<-cbind(apply(cbind((gammas200[,3:4]),gammas200[,5]),2,sd),
            apply(cbind((gammas400[,3:4]),gammas400[,5]),2,sd),
            apply(cbind((gammas600[,3:4]),gammas600[,5]),2,sd),
            apply(cbind((gammas800[,3:4]),gammas800[,5]),2,sd),
            apply(cbind((gammas1000[,3:4]),gammas1000[,5]),2,sd))

#Fazendo os SD dos parametros da parte logistica
dpsg<-cbind(apply(gammas200[,1:2],2,sd),
            apply(gammas400[,1:2],2,sd),
            apply(gammas600[,1:2],2,sd),
            apply(gammas800[,1:2],2,sd),
            apply(gammas1000[,1:2],2,sd))


#Variancia 
varianciag<-cbind(apply((gammas200[,1:2]),2,var),
                  apply((gammas400[,1:2]),2,var),
                  apply((gammas600[,1:2]),2,var),
                  apply((gammas800[,1:2]),2,var),
                  apply((gammas1000[,1:2]),2,var))

viesg<-rbind( logis[1,]+3,logis[2,]-1 )

eqmg<-viesg^2+varianciag  

samp<-c(200,400,600,800,1000)


varg200 <- read.table(paste(wd,"/simulation200/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg400 <- read.table(paste(wd,"/simulation400/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg600 <- read.table(paste(wd,"/simulation600/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg800 <- read.table(paste(wd,"/simulation800/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg1000 <- read.table(paste(wd,"/simulation1000/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))

#Fazendo a média das variancias dos parametros da parte logistica
varsg<-cbind(apply((varg200[,1:2]),2,mean),
             apply((varg400[,1:2]),2,mean),
             apply((varg600[,1:2]),2,mean),
             apply((varg800[,1:2]),2,mean),
             apply((varg1000[,1:2]),2,mean))

varsg1<-cbind(apply(sqrt(varg200[,1:2]),2,mean),
             apply(sqrt(varg400[,1:2]),2,mean),
             apply(sqrt(varg600[,1:2]),2,mean),
             apply(sqrt(varg800[,1:2]),2,mean),
             apply(sqrt(varg1000[,1:2]),2,mean))



#Teorica E(t|x=)
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
mediaw(x=1)
mediaw(x=-1)


##Mean time ZIgamma
po<-NULL
scale.1<-NULL
para<-NULL
mediag<-function(x){
  for (j in 1:5){
    po[j]=exp(logis[1,j]+logis[2,j]*x)/(1+exp(logis[1,j]+logis[2,j]*x))
    scale.1[j]=gammas[2,j]*exp(gammas[3,j]*x)
    para[j]<-gammas[1,j]
  }
  media.g=(1-po)*(para)*(1/scale.1)
  return(rbind(media.g,po))
}

mediag(-1)
mediag(1)

rbind(
mediaw(-1),
mediag(-1),
mediaw(1),
mediag(1))


##Grafico de medias
media2 = data.frame(samp=samp,media.1=c(media1,mediag),grupo=c(rep("ZIWeibull",5),rep("ZIGamma",5)))

ggplot(data=media2, aes(x = samp, y = media.1, group=grupo)) +
  geom_point(show.legend=FALSE, shape = 16) +
  geom_line(aes(linetype = grupo)) +
  geom_hline(yintercept = media) +
  ylim(c(24,27.5)) +
  labs(y = "Mean Time", x = "Sample Size")+
  guides(fill=guide_legend("Regression"))+
  scale_color_discrete(name = "Regression")

#Grafico de po

po2=data.frame(pos=c(po,po1),samp=samp,grupo=c(rep("ZIWeibull",5),rep("ZIGamma",5)))


ggplot(data = po2, aes(x = samp, y = pos, group=grupo)) +
  geom_point(show.legend=FALSE, shape = 16) +
  geom_line(aes(linetype = grupo)) +
  geom_hline(yintercept = p0) +
  ylim(c(min(pos),0.1195)) +
  labs(y = "Mean Time", x = "Sample Size")


#Grafico de betas da logistica

logi<-data.frame(w1=logis[1,],w2=logis[2,],w3=logis[1,]-beta5,w4=logis[2,]-beta6,w5=varianciag[1,],w6=varianciag[2,],samp=samp)

l1 <- ggplot(data = logi, aes(x = samp, y = w1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta5) +
  labs(y = expression(paste("Mean (",hat(beta)[0],")" )), x = "Sample Size")

l2 <- ggplot(data = logi, aes(x = samp, y = w2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta6) +
  labs(y = expression(paste("Mean (",hat(beta)[1],")" )), x = "Sample Size")

l3 <- ggplot(data = logi, aes(x = samp, y = w5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[0], " )")), x = "Sample")

l4 <- ggplot(data = logi, aes(x = samp, y = w4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[1], " )")), x = "Sample")

l5 <- ggplot(data = logi, aes(x = samp, y = w5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[0], " )")), x = "Sample")

l6 <- ggplot(data = logi, aes(x = samp, y = w6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[1], " )")), x = "Sample")


library(gridExtra)
grid.arrange(arrangeGrob(l1,l2,l3,l4,l5,l6, ncol=2))


# --------------------------------------------------------------------------------- #

###Tabelas

##Parte Weibull

##Probabilidade de cobertura

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
final1






##Parte Gamma

estg <- c(-3, 1)

#Criando uma nova matriz que pega apenas os dados gerados para beta_0 e beta_1 (colunas 1 e 2)
gammas200e <- gammas200[,1:2]
gammas400e <- gammas400[,1:2]
gammas600e <- gammas600[,1:2]
gammas800e <- gammas800[,1:2]
gammas1000e <- gammas1000[,1:2]

LimInfg200 <- gammas200e[,]-1.96*sqrt(varianciag[,1])
LimSupg200 <- gammas200e[,]+1.96*sqrt(varianciag[,1])
Rg200 <- nrow(gammas200e)
lim_inf_rg200 <- matrix(NA,Rg200,2)
lim_sup_rg200 <- matrix(NA,Rg200,2)

for (i in 1:2){
  lim_inf_rg200[,i] <- gammas200e[,i] -1.96*sd(gammas200e[,i])
  lim_sup_rg200[,i] <- gammas200e[,i] +1.96*sd(gammas200e[,i])
}

Prob_coberturag200 <- NULL
for (i in 1:2){
  Prob_coberturag200[i] <- (sum(lim_sup_rg200[,i] >= estg[i] & lim_inf_rg200[,i] <= estg[i]))/Rg200
}


LimInfg400 <- gammas400e[,]-1.96*sqrt(varianciag[,1])
LimSupg400 <- gammas400e[,]+1.96*sqrt(varianciag[,1])
Rg400 <- nrow(gammas400e)
lim_inf_rg400 <- matrix(NA,Rg400,2)
lim_sup_rg400 <- matrix(NA,Rg400,2)

for (i in 1:2){
  lim_inf_rg400[,i] <- gammas400e[,i] -1.96*sd(gammas400e[,i])
  lim_sup_rg400[,i] <- gammas400e[,i] +1.96*sd(gammas400e[,i])
}

Prob_coberturag400 <- NULL
for (i in 1:2){
  Prob_coberturag400[i] <- (sum(lim_sup_rg400[,i] >= estg[i] & lim_inf_rg400[,i] <= estg[i]))/Rg400
}


LimInfg600 <- gammas600e[,]-1.96*sqrt(varianciag[,1])
LimSupg600 <- gammas600e[,]+1.96*sqrt(varianciag[,1])
Rg600 <- nrow(gammas600e)
lim_inf_rg600 <- matrix(NA,Rg600,2)
lim_sup_rg600 <- matrix(NA,Rg600,2)

for (i in 1:2){
  lim_inf_rg600[,i] <- gammas600e[,i] -1.96*sd(gammas600e[,i])
  lim_sup_rg600[,i] <- gammas600e[,i] +1.96*sd(gammas600e[,i])
}

Prob_coberturag600 <- NULL
for (i in 1:2){
  Prob_coberturag600[i] <- (sum(lim_sup_rg600[,i] >= estg[i] & lim_inf_rg600[,i] <= estg[i]))/Rg600
}


LimInfg800 <- gammas800e[,]-1.96*sqrt(varianciag[,1])
LimSupg800 <- gammas800e[,]+1.96*sqrt(varianciag[,1])
Rg800 <- nrow(gammas800e)
lim_inf_rg800 <- matrix(NA,Rg800,2)
lim_sup_rg800 <- matrix(NA,Rg800,2)

for (i in 1:2){
  lim_inf_rg800[,i] <- gammas800e[,i] -1.96*sd(gammas800e[,i])
  lim_sup_rg800[,i] <- gammas800e[,i] +1.96*sd(gammas800e[,i])
}

Prob_coberturag800 <- NULL
for (i in 1:2){
  Prob_coberturag800[i] <- (sum(lim_sup_rg800[,i] >= estg[i] & lim_inf_rg800[,i] <= estg[i]))/Rg800
}


LimInfg1000 <- gammas1000e[,]-1.96*sqrt(varianciag[,1])
LimSupg1000 <- gammas1000e[,]+1.96*sqrt(varianciag[,1])
Rg1000 <- nrow(gammas1000e)
lim_inf_rg1000 <- matrix(NA,Rg1000,2)
lim_sup_rg1000 <- matrix(NA,Rg1000,2)

for (i in 1:2){
  lim_inf_rg1000[,i] <- gammas1000e[,i] -1.96*sd(gammas1000e[,i])
  lim_sup_rg1000[,i] <- gammas1000e[,i] +1.96*sd(gammas1000e[,i])
}

Prob_coberturag1000 <- NULL
for (i in 1:2){
  Prob_coberturag1000[i] <- (sum(lim_sup_rg1000[,i] >= estg[i] & lim_inf_rg1000[,i] <= estg[i]))/Rg1000
}

probcobg <- cbind(Prob_coberturag200,
                  Prob_coberturag400,
                  Prob_coberturag600,
                  Prob_coberturag800,
                  Prob_coberturag1000)


#Usaremos logis, que é a matriz que ja faz as medias dos parametros da parte logistica

tabg200<-cbind(logis[,1],cbind(logis[,1]-1.96*sqrt(varianciag[,1]),logis[,1]+1.96*sqrt(varianciag[,1])),dpsg[,1],varsg[,1],viesg[,1],eqmg[,1],probcobg[,1])
row.names(tabg200)<-c(1,2)

tabg400<-cbind(logis[,2],cbind(logis[,2]-1.96*sqrt(varianciag[,2]),logis[,2]+1.96*sqrt(varianciag[,2])),dpsg[,2],varsg[,2],viesg[,2],eqmg[,2],probcobg[,2])
row.names(tabg400)<-c(1,2)

tabg600<-cbind(logis[,3],cbind(logis[,3]-1.96*sqrt(varianciag[,3]),logis[,3]+1.96*sqrt(varianciag[,3])),dpsg[,3],varsg[,3],viesg[,3],eqmg[,3],probcobg[,3])
row.names(tabg600)<-c(1,2)

tabg800<-cbind(logis[,4],cbind(logis[,4]-1.96*sqrt(varianciag[,4]),logis[,4]+1.96*sqrt(varianciag[,4])),dpsg[,4],varsg[,4],viesg[,4],eqmg[,4],probcobg[,4])
row.names(tabg800)<-c(1,2)

tabg1000<-cbind(logis[,5],cbind(logis[,5]-1.96*sqrt(varianciag[,5]),logis[,5]+1.96*sqrt(varianciag[,5])),dpsg[,5],varsg[,5],viesg[,5],eqmg[,5],probcobg[,5])
row.names(tabg1000)<-c(1,2)

final<-rbind(tabg200,tabg400,tabg600,tabg800,tabg1000)
row.names(final)<-rep(c("beta_0", "beta_1"),5)
colnames(final)<-c("media","LCI","UCI","SD of betas","Mean Std Err","vies","eqmg","CP")
final



#Tabela de medias e p0
