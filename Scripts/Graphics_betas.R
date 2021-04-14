#Packages to be used
packages<-c("readxl","here","tidyverse","ggpubr","ggplot2","flexsurv","knitr","glmmsr","plotly","gridExtra","grid","ggridges","ggthemes","summarytools","ggcorrplot","survival")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


wd <- getwd()

###--------------------------------------------------------------------------------------------###
# Graphics for weibull distribution
###--------------------------------------------------------------------------------------------###

# Define files
dir.create(paste(wd,"/Figure_simulation",sep=""))

# Set environment
lime = 10
N<-c(200,400,600,800,1000)
modelos <- c("dados_No_p1", "modelos_No_p1")


##Remember the parameters
beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1
varp=c(beta1,beta2,beta3,beta4,beta5,beta6)

#####Betas ZIWeibull model
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

####Betas ZIGamma model
gammas200 <- read.table(paste(wd,'/','Output/simulation_gamma',N[1],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas400 <- read.table(paste(wd,'/','Output/simulation_gamma',N[2],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas600 <- read.table(paste(wd,'/','Output/simulation_gamma',N[3],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas800 <- read.table(paste(wd,'/','Output/simulation_gamma',N[4],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas1000 <- read.table(paste(wd,'/','Output/simulation_gamma',N[5],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))

#Transform the parameters to natural scale
#V1: beta5=-3
#V2: beta6=1
#V3: beta1=0.5
#V4: beta3=1.5
#V5: beta4=2
#V6: beta2=0.5
#V7: LogLik
#V7: AIC

gammas<-cbind(apply(gammas200[,1:6],2,mean),
              apply(gammas400[,1:6],2,mean),
              apply(gammas600[,1:6],2,mean),
              apply(gammas800[,1:6],2,mean),
              apply(gammas1000[,1:6],2,mean)
)

vies_g<-(gammas-varp)

## SD and Var ZIWeibull model
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

## SD and Var ZIGamma model
dps_g<-cbind(apply(gammas200[,1:6],2,sd),
            apply(gammas400[,1:6],2,sd),
            apply(gammas600[,1:6],2,sd),
            apply(gammas800[,1:6],2,sd),
            apply(gammas1000[,1:6],2,sd))


variancia_g<-cbind(apply(gammas200[,1:6],2,var),
                   apply(gammas400[,1:6],2,var),
                   apply(gammas600[,1:6],2,var),
                   apply(gammas800[,1:6],2,var),
                   apply(gammas1000[,1:6],2,var))

# Standard errors
#Weibull model
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

#Gamma model
varg200 <- read.table(paste(wd,"/Output/simulation_gamma200/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg400 <- read.table(paste(wd,"/Output/simulation_gamma400/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg600 <- read.table(paste(wd,"/Output/simulation_gamma600/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg800 <- read.table(paste(wd,"/Output/simulation_gamma800/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg1000 <- read.table(paste(wd,"/Output/simulation_gamma1000/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))

varg200<-cbind.data.frame(varg200[,1:4],(varg200[,5:6])^2)
varg400<-cbind.data.frame(varg400[,1:4],(varg400[,5:6])^2)
varg600<-cbind.data.frame(varg600[,1:4],(varg600[,5:6])^2)
varg800<-cbind.data.frame(varg800[,1:4],(varg800[,5:6])^2)
varg1000<-cbind.data.frame(varg1000[,1:4],(varg1000[,5:6])^2)

var_g<-cbind(apply((varg200),2,mean),
              apply((varg400),2,mean),
              apply((varg600),2,mean),
              apply((varg800),2,mean),
              apply((varg1000),2,mean)
            )

eqm_g<-(vies_g^2+variancia_g)

#Plots
vies1<-t(vies)
vies1_g = t(vies_g)

vies2 = data.frame(rbind(vies1,vies1_g))
vies2$Model<-c(rep("ZIWeibull",5),rep("ZIGgamma",5))
attach(vies2)
samp<-rep(c(200,400,600,800,1000),2)

a <- ggplot(data = vies2, aes(x = samp, y = V1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[1], " )")), x = "Sample Size")

b <- ggplot(data = vies2, aes(x = samp, y = V2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[2], " )")), x = "Sample Size")

c <- ggplot(data = vies2, aes(x = samp, y = V3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[3], " )")), x = "Sample Size")

d <- ggplot(data = vies2, aes(x = samp, y = V4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[4], " )")), x = "Sample Size")

e <- ggplot(data = vies2, aes(x = samp, y = V5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[5], " )")), x = "Sample Size")

f <- ggplot(data = vies2, aes(x = samp, y = V6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("BIAS ( ", hat(beta)[6], " )")), x = "Sample Size")

tiff(file=paste(wd,"/Figure_simulation","/bias.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
gridExtra::grid.arrange(a,b,c,d,e,f,nrow=3)
dev.off()

rownames(vars) = c("Z1","Z2","Z3","Z4","Z5","Z6") 
rownames(var_g) = c("Z1","Z2","Z3","Z4","Z5","Z6") 
vars1 = t(vars)
vars1_g=t(var_g)

vars2 = data.frame(rbind(vars1,vars1_g))
vars2$Model<-c(rep("ZIWeibull",5),rep("ZIgamma",5))
attach(vars2)
samp<-rep(c(200,400,600,800,1000),2)

vars2 = data.frame(vars2)
attach(vars2)

a2 <- ggplot(data = vars2, aes(x = samp, y = Z1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[1], " )")), x = "Sample Size")+
  theme(legend.position="none")

b2 <- ggplot(data = vars2, aes(x = samp, y = Z2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[2], " )")), x = "Sample Size")+
  theme(legend.position="none")

c2 <- ggplot(data = vars2, aes(x = samp, y = Z3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[3], " )")), x = "Sample Size")+
  theme(legend.position="none")

d2 <- ggplot(data = vars2, aes(x = samp, y = Z4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[4], " )")), x = "Sample Size")+
  theme(legend.position="none")

e2 <- ggplot(data = vars2, aes(x = samp, y = Z5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[5], " )")), x = "Sample Size")+
  theme(legend.position="none")

f2 <- ggplot(data = vars2, aes(x = samp, y = Z6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  labs(y = expression(paste("VAR ( ", hat(beta)[6], " )")), x = "Sample Size")+
  theme(legend.position="none")

tiff(file=paste(wd,"/Figure_simulation","/var.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
ggpubr::ggarrange(a2,b2,c2,d2,e2,f2,nrow=3,ncol=2,common.legend = TRUE, legend="bottom")
dev.off()


#Plots

rownames(betas) = c("W1","W2","W3","W4","W5","W6")
rownames(gammas) = c("W1","W2","W3","W4","W5","W6")
betas1 = t(betas)
betas1_g = t(gammas)
betas2 = data.frame(rbind(betas1,betas1_g))
betas2$Model<-c(rep("ZIWeibull",5),rep("ZIGamma",5))
attach(betas2)
samp<-rep(c(200,400,600,800,1000),2)

a3 <- ggplot(data = betas2, aes(x = samp, y = W1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype = Model)) +
  geom_hline(yintercept = beta1) +
  labs(y = expression(paste("Mean (",hat(beta)[1],")" )), x = "Sample Size")+
  theme(legend.position="none")

b3 <- ggplot(data = betas2, aes(x = samp, y = W2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype =  Model)) +
  geom_hline(yintercept = beta2) +
  labs(y = expression(paste("Mean (",hat(beta)[2],")" )), x = "Sample Size")+
  theme(legend.position="none")

c3 <- ggplot(data = betas2, aes(x = samp, y = W3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype =  Model)) +
  geom_hline(yintercept = beta3) +
  labs(y = expression(paste("Mean (",hat(beta)[3],")" )), x = "Sample Size")+
  theme(legend.position="none")

d3 <- ggplot(data = betas2, aes(x = samp, y = W4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype =  Model)) +
  geom_hline(yintercept = beta4) +
  labs(y = expression(paste("Mean (",hat(beta)[4],")" )), x = "Sample Size")+
  theme(legend.position="none")

e3 <- ggplot(data = betas2, aes(x = samp, y = W5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype =  Model)) +
  geom_hline(yintercept = beta5) +
  labs(y = expression(paste("Mean (",hat(beta)[5],")" )), x = "Sample Size")+
  theme(legend.position="none")

f3 <- ggplot(data = betas2, aes(x = samp, y = W6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(aes(linetype =  Model)) +
  geom_hline(yintercept = beta6) +
  labs(y = expression(paste("Mean (",hat(beta)[6],")" )), x = "Sample Size")+
  theme(legend.position="none")

tiff(file=paste(wd,"/Figure_simulation","/beta.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
ggpubr::ggarrange(a3,b3,c3,d3,e3,f3,nrow=3,ncol=2,common.legend = TRUE, legend="bottom")

#ggarrange(arrangeGrob(a3,b3, ncol=2), 
#             arrangeGrob(c3,d3, ncol=2),
#             arrangeGrob(e3,f3, ncol=2), ncol=1)

dev.off()




###--------------------------------------------------------------------------------------------###
# Tables
###--------------------------------------------------------------------------------------------###

# Define files
dir.create(paste(wd,"/Tables_simulation",sep=""))

## Theoretic E(t|x=) for Weibull
media<-function(x){
  alfa=exp(beta1+beta2*x)
  teta=exp(beta3+beta4*x)
  p0=exp(beta5+beta6*x)/(1+exp(beta5+beta6*x))
  
  media.1=(1-p0)*teta*gamma(1+(1/alfa))
  return(rbind(media.1,p0))
}

media(-1);media(1)

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

means_w<-rbind(
  mediaw(x=1),
  mediaw(x=-1))

colnames(means_w)[1:5]<-c("200","400","600","800","1000")

row.names(means_w)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(means_w,here("Tables_simulation","means_weibull.txt"),sep=";",row.names = T)




## Theoretic E(t|x=) for gamma
media1<-function(x){
  alfa=exp(beta1+beta2*x)
  teta=exp(beta3+beta4*x)
  p0=exp(beta5+beta6*x)/(1+exp(beta5+beta6*x))
  
  media.1=(1-p0)*teta*alfa
  return(rbind(media.1,p0))
}

media1(-1);media1(1)



#Mean time ZIGamma

mediag<-function(x){
  alfa1=numeric()
  teta1=numeric()
  po1=numeric()
  
  for (j in 1:5){
    alfa1[j]=exp(gammas[1,j]+gammas[2,j]*x) #k in the article
    teta1[j]=exp(gammas[3,j]+gammas[4,j]*x) #lambda in the article
    po1[j]=exp(gammas[5,j]+gammas[6,j]*x)/(1+exp(gammas[5,j]+gammas[6,j]*x))
  }

  media1=(1-po1)*teta1*(alfa1)
  return (rbind(media1,po1))
}


means_g<-rbind(
  mediag(x=1),
  mediag(x=-1))

colnames(means_g)[1:5]<-c("200","400","600","800","1000")

row.names(means_g)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(means_g,here("Tables_simulation","means_gamma.txt"),sep=";",row.names = T)



##########
###Tables#
##########

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


tab200<-cbind(betas[,1],cbind(betas[,1]-1.96*sqrt(variancia[,1]),betas[,1]+1.96*sqrt(variancia[,1])),dps[,1],sqrt(var_g[,1]),vies_g[,1],eqm_g[,1],probcob[,1])
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

write.table(final1,here("Tables_simulation","final_weibull.txt"),sep=";",row.names = F)






##Gamma part

##Cover probability

gammas200e <- gammas200[,1:6]
gammas400e <- gammas400[,1:6]
gammas600e <- gammas600[,1:6]
gammas800e <- gammas800[,1:6]
gammas1000e <- gammas1000[,1:6]

LimInf200 <- gammas200e[,]-1.96*sqrt(variancia_g[,1])
LimSup200 <- gammas200e[,]+1.96*sqrt(variancia_g[,1])
R200 <- nrow(gammas200e)
lim_inf_r200 <- matrix(NA,R200,6)
lim_sup_r200 <- matrix(NA,R200,6)

for (i in 1:6){
  lim_inf_r200[,i] <- gammas200e[,i] -1.96*sd(gammas200e[,i])
  lim_sup_r200[,i] <- gammas200e[,i] +1.96*sd(gammas200e[,i])
}
Prob_cobertura200 <- NULL
for (i in 1:6){
  Prob_cobertura200[i] <- (sum(lim_sup_r200[,i] >= varp[i] & lim_inf_r200[,i] <= varp[i]))/R200
}


LimInf400 <- gammas400e[,]-1.96*sqrt(variancia_g[,2])
LimSup400 <- gammas400e[,]+1.96*sqrt(variancia_g[,2])
R400 <- nrow(gammas400e)
lim_inf_r400 <- matrix(NA,R400,6)
lim_sup_r400 <- matrix(NA,R400,6)

for (i in 1:6){
  lim_inf_r400[,i] <- gammas400e[,i] -1.96*sd(gammas400e[,i])
  lim_sup_r400[,i] <- gammas400e[,i] +1.96*sd(gammas400e[,i])
}
Prob_cobertura400 <- NULL
for (i in 1:6){
  Prob_cobertura400[i] <- (sum(lim_sup_r400[,i] >= varp[i] & lim_inf_r400[,i] <= varp[i]))/R400
}


LimInf600 <- gammas600e[,]-1.96*sqrt(variancia_g[,3])
LimSup600 <- gammas600e[,]+1.96*sqrt(variancia_g[,3])
R600 <- nrow(gammas600e)
lim_inf_r600 <- matrix(NA,R600,6)
lim_sup_r600 <- matrix(NA,R600,6)

for (i in 1:6){
  lim_inf_r600[,i] <- gammas600e[,i] -1.96*sd(gammas600e[,i])
  lim_sup_r600[,i] <- gammas600e[,i] +1.96*sd(gammas600e[,i])
}
Prob_cobertura600 <- NULL
for (i in 1:6){
  Prob_cobertura600[i] <- (sum(lim_sup_r600[,i] >= varp[i] & lim_inf_r600[,i] <= varp[i]))/R600
}


LimInf800 <- gammas800e[,]-1.96*sqrt(variancia_g[,4])
LimSup800 <- gammas800e[,]+1.96*sqrt(variancia_g[,4])
R800 <- nrow(gammas800e)
lim_inf_r800 <- matrix(NA,R800,6)
lim_sup_r800 <- matrix(NA,R800,6)

for (i in 1:6){
  lim_inf_r800[,i] <- gammas800e[,i] -1.96*sd(gammas800e[,i])
  lim_sup_r800[,i] <- gammas800e[,i] +1.96*sd(gammas800e[,i])
}
Prob_cobertura800 <- NULL
for (i in 1:6){
  Prob_cobertura800[i] <- (sum(lim_sup_r800[,i] >= varp[i] & lim_inf_r800[,i] <= varp[i]))/R800
}

LimInf1000 <- gammas1000e[,]-1.96*sqrt(variancia_g[,5])
LimSup1000 <- gammas1000e[,]+1.96*sqrt(variancia_g[,5])
R1000 <- nrow(gammas1000e)
lim_inf_r1000 <- matrix(NA,R1000,6)
lim_sup_r1000 <- matrix(NA,R1000,6)

for (i in 1:6){
  lim_inf_r1000[,i] <- gammas1000e[,i] -1.96*sd(gammas1000e[,i])
  lim_sup_r1000[,i] <- gammas1000e[,i] +1.96*sd(gammas1000e[,i])
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


tab200<-cbind(gammas[,1],cbind(gammas[,1]-1.96*sqrt(variancia_g[,1]),gammas[,1]+1.96*sqrt(variancia_g[,1])),dps_g[,1],sqrt(var_g[,1]),vies_g[,1],eqm_g[,1],probcob[,1])
row.names(tab200)<-c(1,2,3,4,5,6)

tab400<-cbind(gammas[,2],cbind(gammas[,2]-1.96*sqrt(variancia_g[,2]),gammas[,2]+1.96*sqrt(variancia_g[,2])),dps_g[,2],sqrt(var_g[,2]),vies_g[,2],eqm_g[,2],probcob[,2])
row.names(tab400)<-c(1,2,3,4,5,6)


tab600<-cbind(gammas[,3],cbind(gammas[,3]-1.96*sqrt(variancia_g[,3]),gammas[,3]+1.96*sqrt(variancia_g[,3])),dps_g[,3],sqrt(var_g[,3]),vies_g[,3],eqm_g[,3],probcob[,3])
row.names(tab600)<-c(1,2,3,4,5,6)

tab800<-cbind(gammas[,4],cbind(gammas[,4]-1.96*sqrt(variancia_g[,4]),gammas[,4]+1.96*sqrt(variancia_g[,4])),dps_g[,4],sqrt(var_g[,4]),vies_g[,4],eqm_g[,4],probcob[,4])
row.names(tab800)<-c(1,2,3,4,5,6)


tab1000<-cbind(gammas[,5],cbind(gammas[,5]-1.96*sqrt(variancia_g[,5]),gammas[,5]+1.96*sqrt(variancia_g[,5])),dps_g[,5],sqrt(var_g[,5]),vies_g[,5],eqm_g[,5],probcob[,5])
row.names(tab1000)<-c(1,2,3,4,5,6)



final2<-rbind(tab200,tab400,tab600,tab800,tab1000)
row.names(final2)<-rep(c(1:6),5)
colnames(final2)<-c("media","LCI","UCI","SD of betas","Mean Std Err","vies","eqmg","CP")

write.table(final2,here("Tables_simulation","final_gamma.txt"),sep=";",row.names = F)
