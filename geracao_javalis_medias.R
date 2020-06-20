set.seed(13)
library(knitr)
library(ggplot2)


#Sample size
n <- 400

#Criando sexo
z<- rep(1,n/2) 
w<-rep(0,n/2)
x1<-  c(z,w)
  
#Criando idades em meses
age_fem<-rexp(n/2,1/12)
age_mac<-rexp(n/2,1/15)
ages<-c(age_fem,age_mac)
  
# Criando distancias, zero e censura. 
zer<- (exp(-1.5 - 0.02*(ages) + 0.2*x1))/(1+exp(-1.5 - 0.02*(ages) + 0.2*x1))  #Probabilidade de ser zero


alpha.t <- 1 #faz com que fique exponencial
lambda.t <- 0.1

#fem
gamma.t <- lambda.t*exp(-0.02*(age_fem)+ 0.5*z)
u <- runif(n)
y_fem_auxi <- NULL

for(i in 1:(n/2)){
  y_fem_auxi[i] <- (-log(u[i])/gamma.t[i])^(1/alpha.t) 
}

y_fem <- y_fem_auxi*(1-rbinom(200,1,zer[1:200]))  #Distancias para femeas

summary(y_fem)
hist(y_fem)

#mac
gamma.t <- lambda.t*exp(-0.02*(age_mac) + 0.5*w)
u <- runif(n)
y_mac_auxi <- NULL

for(i in 1:(n/2)){
  y_mac_auxi[i] <- (-log(u[i])/gamma.t[i])^(1/alpha.t) 
}

y_mac <- y_mac_auxi*(1-rbinom(200,1,zer[201:400])) #Distancias para machos

summary(y_mac)
hist(y_mac)

y<-c(y_fem, y_mac ) 
   
c <- rep(30,n)
   
t <- NULL
delta <- NULL
Z <- NULL
for (i in 1:n){ 
    t[i] = min(y[i],c[i]) #tempo observado (minimo entre os tempos de falha e censura)
  
      if (t[i] < c[i]){
        delta[i]=1          #delta: indicadora de falha e censura
      } else {delta[i]=0}
  
      if (t[i] == 0){ #Z: indicadora de zero
        Z[i]=1
      }else{Z[i]=0}
  }
   
#Criando data frame
data <- as.data.frame(cbind(t,Z,delta, ages,x1))
names(data) <- c("dist","zero","delta","age", "sex")  

#percentual de eventos
sum(data$delta)/n
#[1] 0.9725

#percentual de zeros
sum(data$zero)/n
#[1] 0.2375

require(survival)
plot(Surv(data$dist,delta))


##################
#Resumo dos dados#
##################

tabela<-cbind(
  #Dist?ncia media geral e por sexo
  rbind(
    mean(data$dist[data$sex==0]), # geral para machos
    mean(data$dist[data$sex==0&data$age<15]), #Machos <16
    mean(data$dist[data$sex==0&data$age>15]), #Machos >16
    mean(data$dist[data$sex==1]), # geral para femeas
    mean(data$dist[data$sex==1&data$age<12]), #Femeas <16
    mean(data$dist[data$sex==1&data$age>12]) #Femeas >16
  ),
  
  # %zeros media geral e por sexo
  rbind(
    mean(data$zero[data$sex==0]), # geral para machos
    mean(data$zero[data$sex==0&data$age<15]), #Machos <16
    mean(data$zero[data$sex==0&data$age>15]), #Machos >16
    mean(data$zero[data$sex==1]), # geral para femeas
    mean(data$zero[data$sex==1&data$age<12]), #Femeas <16
    mean(data$zero[data$sex==1&data$age>12]) #Femeas >16
  ),
  
  
  #% censura media geral e por sexo
  rbind(
    mean(data$delta[data$sex==0]), # geral para machos
    mean(data$delta[data$sex==0&data$age<15]), #Machos <16
    mean(data$delta[data$sex==0&data$age>15]), #Machos >16
    mean(data$delta[data$sex==1]), # geral para femeas
    mean(data$delta[data$sex==1&data$age<12]), #Femeas <16
    mean(data$delta[data$sex==1&data$age>12]) #Femeas >16
  )
  
)

row.names(tabela)<-c("Machos","MAchos <16","MAchos >16","Femeas","Femeas<16","Femeas>16")
colnames(tabela)<-c("Dist media","%zeros","%falhas")
kable(tabela)

par(mfrow=c(1,2))
hist(data$dist[data$sex==1],main="Histogram of distance for females",xlab="Distance")
hist(data$dist[data$sex==0],main="Histogram of distance for males",xlab="Distance")

data$Sex<-c(rep("Female",200),rep("Male",200))
p<-ggplot(data, aes(x=dist))+
  geom_histogram(color="black", fill="white")+
  facet_grid(Sex ~ .)+
  labs(x ="Distance (Km)", y = "Frequency")
p

##################
# Modelo Weibull #
##################


GPE=function(varp)
{
  a       =   exp(varp[1])*exp(x1*varp[2])*exp(x2*varp[3])
  lam     =   exp(varp[4])*exp(x1*varp[5])*exp(x2*varp[6])
  #gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  gamm_zer = (exp(varp[7])*exp(x1*varp[8])*exp(x2*varp[9]))/(1+exp(varp[7])*exp(x1*varp[8])*exp(x2*varp[9]))
  #gamm_inf = (exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  f0   = gamm_zer
  fw = (exp(log(a)-log(lam)))*((exp(log(time)-log(lam)))^(a-1))*(exp(-(exp(log(time)-log(lam)))^a))
  Sw = (exp(-(exp(log(time)-log(lam)))^a))
  #fpop = (1-gamm_zer-gamm_inf)*fw
  #Spop = gamm_inf + (1-gamm_zer-gamm_inf)*Sw
  fpop = (1-gamm_zer)*fw
  Spop = (1-gamm_zer)*Sw
  f    = (fpop^delta)*(Spop^(1-delta))
  g   <- ifelse(time==0, f0, f)
  adFunc   = sum(log(g))
  return(adFunc)
}


# chute inicial 
varp=rep(0,9)

#dados observados
x1    <- (data$age)
x2    <- data$sex
time  <- data$dist
delta <- data$delta

fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))

estc = try(fit$par)
Hc=try(fit$hessian); varic=try(-solve(Hc))


cbind(estc,sqrt(diag(varic)),estc-1.96*sqrt(diag(varic)),estc+1.96*sqrt(diag(varic)))



#Modelo Gamma
#############

library(flexsurv)


summary(glm(data$zero ~ (data$age)+data$sex,family=binomial(link='logit'))->a)

b<-flexsurvreg(Surv(data$dist, data$delta,type='right')~(data$age)+data$sex,dist="gamma",subset = data$zer==0)


##Media modelo Weibull

#Distancia esperada 
#x=age; y=sex{1=femea, 0=macho}

mw1<-function(x,y){
  p0<-(exp(estc[7]+estc[8]*x+estc[9]*y)/(1+exp(estc[7]+estc[8]*x+estc[9]*y)))
  alpha<-(exp(estc[1]+estc[2]*x+estc[3]*y))
  theta<-exp(estc[4]+estc[5]*x+estc[6]*y)
  mw1 <- (1-p0)*theta*gamma(1+1/alpha)
  mw1
}
#femeas
mw1(12,1)
#machos
mw1(15,0)


###Media modelo gamma
#Distancia esperada 
#x=age; y=sex{1=femea, 0=macho}

mg1<-function(x,y){
  scale.1<-( ((b$res[2]))*exp(b$res[3]*x+b$res[4]*y))
  p=exp(a$coefficients[1]+a$coefficients[2]*x+a$coefficients[3]*y)/(1+(exp(a$coefficients[1]+a$coefficients[2]*x+a$coefficients[3]*y)))
  mg1<-(1-p)*b$res[1]*(1/scale.1)
  mg1
}
#femeas
mg1(12,1)
#machos
mg1(15,0)



##Plots
dom<-seq(0,20,0.1)

media<-c(mw1(dom,1),
         mw1(dom,0),
         mg1(dom,1),
         mg1(dom,0))


#Mean for the mean age
media2<-c(mw1(mean(data$age[data$sex==1]),1),
          mw1(mean(data$age[data$sex==0]),0),
          mg1(mean(data$age[data$sex==1]),1),
          mg1(mean(data$age[data$sex==0]),0))



Sex<-c(rep("Female",length(dom)),rep("Male",length(dom)),rep("Female",length(dom)),rep("Male",length(dom)))

Regression<-c(rep("ZIWeibull",length(dom)*2),rep("ZIGamma",length(dom)*2))

media2 = data.frame(Scaled_Age=c(dom,dom,dom,dom),media.1=media,Sex=Sex,grupo=Regression)

ggplot(data=media2, aes(x = Scaled_Age, y = media.1, group=Sex)) +
#  geom_point(show.legend=FALSE, shape = 10) +
  geom_line(aes(linetype = Sex)) +
  ylim(c(0,max(mw1(dom,0)))) +
  labs(y = "Mean Distance Km", x = "Age in month")+
  guides(fill=guide_legend("Regression"))+
  scale_color_discrete(name = "Regression")+
  facet_grid(grupo~.)


#testes de hipoteses weibull
betas_weibull=estc
SE_weibull=sqrt(diag(varic))
wald_weibull=betas_weibull^2/SE_weibull^2
pchisq(wald_weibull,1,lower.tail = F)
Lic=betas_weibull-1.96*SE_weibull
Uic=betas_weibull+1.96*SE_weibull

#testes de hipoteses gamma
betas_gamma=b$res[,1]
SE_gamma=b$res[,4]
wald_gamma = betas_gamma^2/SE_gamma^2
pchisq(wald_gamma,1,lower.tail = F)


##Criterios de informacao

AIC_gamma=(2*7)-2*(logLik(b)+logLik(a))
AIC_weibull=(2*9)-2*fit$value
AICc_gamma=AIC_gamma-(2*7*(7+1)/(400-7-1))
AICc_weibull=AIC_weibull-(2*9*(9+1)/(400-9-1))
BIC_gamma=(log(400)*7)-2*(logLik(b)+logLik(a))
BIC_weibull=(log(400)*9)-2*fit$value

AIC_gamma
AIC_weibull
AICc_gamma
AICc_weibull
BIC_gamma
BIC_weibull


##Modelo linear

hist(data$dist[data$zero==0 & data$delta==1])

tra<-data$dist[data$zero==0 & data$delta==1]

data2<-subset(data,data$zero==0 & data$delta==1)


summary(lm(log(dist)~sex+age,data=data2)->mod)
hist(mod$residuals)


plot(mod)

exp(mod$coefficients[1]+mod$coefficients[2]*0+mod$coefficients[3]*15)

exp(mod$coefficients[1]+mod$coefficients[2]*1+mod$coefficients[3]*12)

