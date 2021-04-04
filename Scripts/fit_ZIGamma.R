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
    
    c<-flexsurvreg(Surv(t, d,type='right')~x,
                   anc = list(shape = ~ x),data=dados,dist="gamma",subset=dados$z==0)
    
    
    
    
    varic=try(solve(c$opt$hessian))
    
    wd.sim<-paste(wd,"/Output/simulation",N[j],"/modelos_No_p1/", 'gamma', sep="")  
    param_Gamma<- file.path(wd.sim, "param_Gamma.txt") 
    Erro_Gamma <- file.path(wd.sim, "Erro_Gamma.txt")
    
    cat(c(a$coefficients,c$res[1],c$res[2],c$res[3],c$res[4],c$loglik,c$AIC),file = param_Gamma,append = TRUE, "\n") 
    cat(c(b$coefficients[1,2],b$coefficients[2,2],varic[1,1],varic[2,2],varic[3,3],varic[4,4]),file = Erro_Gamma,append = TRUE, "\n") 
  }
  
}


dir.create(paste(wd,"/Figure_gamma",sep=""))

#Gamma model

##Graphics
gammas200 <- read.table(paste(wd,'/','Output/simulation',N[1],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas400 <- read.table(paste(wd,'/','Output/simulation',N[2],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas600 <- read.table(paste(wd,'/','Output/simulation',N[3],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas800 <- read.table(paste(wd,'/','Output/simulation',N[4],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas1000 <- read.table(paste(wd,'/','Output/simulation',N[5],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))

#exp means
gammas<-cbind(apply(gammas200[,3:6],2,mean),
              apply(gammas400[,3:6],2,mean),
              apply(gammas600[,3:6],2,mean),
              apply(gammas800[,3:6],2,mean),
              apply(gammas1000[,3:6],2,mean)
        )

#means logistic part
logis<-cbind(apply((gammas200[,1:2]),2,mean),
             apply((gammas400[,1:2]),2,mean),
             apply((gammas600[,1:2]),2,mean),
             apply((gammas800[,1:2]),2,mean),
             apply((gammas1000[,1:2]),2,mean)
      )

#SD of the parameters (exp)
gdps<-cbind(apply(gammas200[,3:6],2,sd),
            apply(gammas400[,3:6],2,sd),
            apply(gammas600[,3:6],2,sd),
            apply(gammas800[,3:6],2,sd),
            apply(gammas1000[,3:6],2,sd)
      )

#SD of parameters logistic part
dpsg<-cbind(apply(gammas200[,1:2],2,sd),
            apply(gammas400[,1:2],2,sd),
            apply(gammas600[,1:2],2,sd),
            apply(gammas800[,1:2],2,sd),
            apply(gammas1000[,1:2],2,sd)
      )


#Variance
varianciag<-cbind(apply((gammas200[,1:6]),2,var),
                  apply((gammas400[,1:6]),2,var),
                  apply((gammas600[,1:6]),2,var),
                  apply((gammas800[,1:6]),2,var),
                  apply((gammas1000[,1:6]),2,var)
            )

viesg<-rbind( logis[1,]+3,logis[2,]-1 )


eqmg<-viesg^2+varianciag  

samp<-c(200,400,600,800,1000)


varg200 <- read.table(paste(wd,"/Output/simulation200/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg400 <- read.table(paste(wd,"/Output/simulation400/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg600 <- read.table(paste(wd,"/Output/simulation600/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg800 <- read.table(paste(wd,"/Output/simulation800/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg1000 <- read.table(paste(wd,"/Output/simulation1000/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))

#Mean variance logistic part
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



#Plots of betas parameters from gamma and logit models

logi<-data.frame(w1=logis[1,],w2=logis[2,],w3=gammas[1,],w4=gammas[2,],w5=gammas[3,],
                 w6=gammas[4,],samp=samp)

l1 <- ggplot(data = logi, aes(x = samp, y = w3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("Mean ( ", hat(beta)[1], " )")), x = "Sample Size")

l2 <- ggplot(data = logi, aes(x = samp, y = w4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("Mean ( ", hat(beta)[2], " )","Rate")), x = "Sample Size")

l3 <- ggplot(data = logi, aes(x = samp, y = w5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("Mean ( ", hat(beta)[3], " )")), x = "Sample Size")

l4 <- ggplot(data = logi, aes(x = samp, y = w6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("Mean ( ", hat(beta)[4], " )")), x = "Sample Size")

l5 <- ggplot(data = logi, aes(x = samp, y = w1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta5) +
  labs(y = expression(paste("Mean (",hat(beta)[5],")" )), x = "Sample Size")

l6 <- ggplot(data = logi, aes(x = samp, y = w2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  geom_hline(yintercept = beta6) +
  labs(y = expression(paste("Mean (",hat(beta)[6],")" )), x = "Sample Size")


png(file=paste(wd,"/Figure_gamma","/mean_gamma.png",sep=""), height = 4, width = 6, units = 'in', res=300)
grid.arrange(arrangeGrob(l1,l2,l3,l4,l5,l6, ncol=2))
dev.off()


#Plots of variances from gamma and logit models

varg<-varianciag

logi1<-data.frame(w1=varg[1,],w2=varg[2,],w3=varg[3,],w4=varg[4,],w5=varg[5,],
                 w6=varg[6,],samp=samp)

v1 <- ggplot(data = logi1, aes(x = samp, y = w3)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[1], " )")), x = "Sample Size")

v2 <- ggplot(data = logi1, aes(x = samp, y = w4)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[2], " )")), x = "Sample Size")

v3 <- ggplot(data = logi1, aes(x = samp, y = w5)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[3], " )")), x = "Sample Size")

v4 <- ggplot(data = logi1, aes(x = samp, y = w6)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR ( ", hat(beta)[4], " )")), x = "Sample Size")

v5 <- ggplot(data = logi1, aes(x = samp, y = w1)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR (",hat(beta)[5],")" )), x = "Sample Size")

v6 <- ggplot(data = logi1, aes(x = samp, y = w2)) +
  geom_point(show.legend = FALSE, shape = 16) +
  geom_line(linetype = 2) +
  labs(y = expression(paste("VAR (",hat(beta)[6],")" )), x = "Sample Size")


png(file=paste(wd,"/Figure_gamma","/var_gamma.png",sep=""), height = 4, width = 6, units = 'in', res=300)
grid.arrange(arrangeGrob(v1,v2,v3,v4,v5,v6, ncol=2))
dev.off()

##Tables
dir.create(paste(wd,"/Tables_gamma",sep=""))


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

##Mean time ZIgamma
po<-NULL
scale.1<-NULL
para<-NULL
mediag<-function(x){
  for (j in 1:5){
    po[j]=exp(logis[1,j]+logis[2,j]*x)/(1+exp(logis[1,j]+logis[2,j]*x))
    scale.1[j]=gammas[2,j]*exp(gammas[3,j]*x)
    para[j]<-gammas[1,j]*exp(gammas[4,j]*x)
  }
  media.g=(1-po)*(para)*(1/scale.1)
  return(rbind(media.g,po))
}

means<-rbind(
mediag(-1),
mediag(1)
)


colnames(means)[1:5]<-c("200","400","600","800","1000")

row.names(means)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(round(means,3),here("Tables_gamma","means.txt"),sep=";",row.names = T)


##Gamma part

estg <- c(-3, 1)

#MAtrix only for beta_0 e beta_1 (col 1 e 2)
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


#We will use the logis, which is the matrix with the mean of the parameters from the logistic part

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

kable(round(final,3))

write.table(round(final,3),here("Tables_gamma","final.txt"),sep=";",row.names = F)


