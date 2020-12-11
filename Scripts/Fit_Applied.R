
##################
# Weibull model  #
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


# initial value
varp=rep(0,9)

#observed data
x1    <- (data$age)
x2    <- data$sex
time  <- data$dist
delta <- data$delta

fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))

estc = try(fit$par)
Hc=try(fit$hessian); varic=try(-solve(Hc))


weibull_parameters<-cbind(estc,sqrt(diag(varic)),estc-1.96*sqrt(diag(varic)),estc+1.96*sqrt(diag(varic)))

colnames(weibull_parameters)<-c("Parameter","Standar error","2.5%","97.5%")

write.table(weibull_parameters,here("Tables_applied","weibull.txt"),sep=";",row.names = F)

##############
#Gamma model #
##############

summary(glm(data$zero ~ (data$age)+data$sex,family=binomial(link='logit'))->a)

b<-flexsurvreg(Surv(data$dist, data$delta,type='right')~(data$age)+data$sex,dist="gamma",subset = data$zer==0)


##Mean weibull model

#Expected distance
#x=age; y=sex{1=female, 0=male}

mw1<-function(x,y){
  p0<-(exp(estc[7]+estc[8]*x+estc[9]*y)/(1+exp(estc[7]+estc[8]*x+estc[9]*y)))
  alpha<-(exp(estc[1]+estc[2]*x+estc[3]*y))
  theta<-exp(estc[4]+estc[5]*x+estc[6]*y)
  mw1 <- (1-p0)*theta*gamma(1+1/alpha)
  mw1
}
#female
mw1(12,1)
#male
mw1(15,0)


###Mean distance gamma
#Expected distance 
#x=age; y=sex{1=female, 0=male}

mg1<-function(x,y){
  scale.1<-( ((b$res[2]))*exp(b$res[3]*x+b$res[4]*y))
  p=exp(a$coefficients[1]+a$coefficients[2]*x+a$coefficients[3]*y)/(1+(exp(a$coefficients[1]+a$coefficients[2]*x+a$coefficients[3]*y)))
  mg1<-(1-p)*b$res[1]*(1/scale.1)
  mg1
}
#female
mg1(12,1)
#male
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


tiff(file=paste(wd,"/Figure_applied","/mean_distance.tiff",sep=""), height = 4, width = 6, units = 'in', res=300)
par(xpd=NA)
ggplot(data=media2, aes(x = Scaled_Age, y = media.1, group=Sex)) +
  #  geom_point(show.legend=FALSE, shape = 10) +
  geom_line(aes(linetype = Sex)) +
  ylim(c(0,max(mw1(dom,0)))) +
  labs(y = "Mean Distance Km", x = "Age in month")+
  guides(fill=guide_legend("Regression"))+
  scale_color_discrete(name = "Regression")+
  facet_grid(grupo~.)
dev.off()

#Hypothesis test weibull
betas_weibull=estc
SE_weibull=sqrt(diag(varic))
wald_weibull=betas_weibull^2/SE_weibull^2
pchisq(wald_weibull,1,lower.tail = F)
Lic=betas_weibull-1.96*SE_weibull
Uic=betas_weibull+1.96*SE_weibull

#Hypothesis test gamma
betas_gamma=b$res[,1]
SE_gamma=b$res[,4]
wald_gamma = betas_gamma^2/SE_gamma^2
pchisq(wald_gamma,1,lower.tail = F)


##Information criteria

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


##Linear model

#hist(data$dist[data$zero==0 & data$delta==1])

#tra<-data$dist[data$zero==0 & data$delta==1]

#data2<-subset(data,data$zero==0 & data$delta==1)


#summary(lm(log(dist)~sex+age,data=data2)->mod)
#hist(mod$residuals)


#plot(mod)

#exp(mod$coefficients[1]+mod$coefficients[2]*0+mod$coefficients[3]*15)

#exp(mod$coefficients[1]+mod$coefficients[2]*1+mod$coefficients[3]*12)

