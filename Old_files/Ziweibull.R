##################
# Weibull model  #
##################



ZIweibull<-function(t,d,data,co=c(),dist){

    attach(data)

  if(dist=="weibull"){

    if(length(co)==0){

          GPE=function(varp){
          a       =   exp(varp[1])
          lam     =   exp(varp[2])

          gamm_zer = (exp(varp[3]))/(1+exp(varp[3]))

          f0   = gamm_zer
          fw = (exp(log(a)-log(lam)))*((exp(log(t)-log(lam)))^(a-1))*(exp(-(exp(log(t)-log(lam)))^a))
          Sw = (exp(-(exp(log(t)-log(lam)))^a))

          fpop = (1-gamm_zer)*fw
          Spop = (1-gamm_zer)*Sw
          f    = (fpop^d)*(Spop^(1-d))
          g   <- ifelse(t==0, f0, f)
          adFunc   = sum(log(g))
          return(adFunc)
        }


# initial value
varp=rep(0,3)

fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))

estc = try(fit$par)
Hc=try(fit$hessian); varic=try(-solve(Hc))



#Hypothesis test weibull
betas_weibull=estc
SE_weibull=sqrt(diag(varic))
wald_weibull=betas_weibull^2/SE_weibull^2


p_val<-numeric()
for (i in 1:length(fit$par)){
if(pchisq(wald_weibull[i],1,lower.tail = F)<0.001){
  p_val[i]<-"<0.001"
}else p_val[i]<-pchisq(wald_weibull[i],1,lower.tail = F)
}


AIC=(2*length(fit$par))-2*fit$value
AICc=AIC-(2*length(fit$par)*(length(fit$par)+1)/(length(t)-length(fit$par)-1))
BIC=(log(length(t))*length(fit$par))-2*fit$value

main<-cbind.data.frame(estc,sqrt(diag(varic)),p_val)
colnames(main)<-c("Estimates","Std Err","P-value")
rownames(main)<-c("Shape",co,"Scale",co,"Intercept",co)

cal<- (paste("ZIweibull:t ~ "))


out<-list(cal,main,fit$value,length(t),sum(t==0),sum(d==0),length(fit$par)-1,AIC,AICc,BIC,wild_boar,fit$hessian)
names(out)<-c("Call","Coefficients","LogLik","N","Zeros","Censored","df","AIC","AICc","BIC","Data","Hessian")

return(out)


  }else{

      if(length(co)==2){

        attach(data)

        data_w<-cbind.data.frame(data[, c(co)])

        x1<-data_w[,1]
        x2<-data_w[,2]


       GPE=function(varp){
        a       =   exp(varp[1])*exp(x1*varp[2])*exp(x2*varp[3])
        lam     =   exp(varp[4])*exp(x1*varp[5])*exp(x2*varp[6])
        gamm_zer = (exp(varp[7])*exp(x1*varp[8])*exp(x2*varp[9]))/(1+exp(varp[7])*exp(x1*varp[8])*exp(x2*varp[9]))

        f0   = gamm_zer
        fw = (exp(log(a)-log(lam)))*((exp(log(t)-log(lam)))^(a-1))*(exp(-(exp(log(t)-log(lam)))^a))
        Sw = (exp(-(exp(log(t)-log(lam)))^a))


        fpop = (1-gamm_zer)*fw
        Spop = (1-gamm_zer)*Sw
        f    = (fpop^d)*(Spop^(1-d))
        g   <- ifelse(t==0, f0, f)
        adFunc   = sum(log(g))
        return(adFunc)
      }


      # initial value
      varp=rep(0,9)

      fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))

      estc = try(fit$par)
      Hc=try(fit$hessian); varic=try(-solve(Hc))



      #Hypothesis test weibull
      betas_weibull=estc
      SE_weibull=sqrt(diag(varic))
      wald_weibull=betas_weibull^2/SE_weibull^2


      p_val<-numeric()
      for (i in 1:length(fit$par)){
        if(pchisq(wald_weibull[i],1,lower.tail = F)<0.001){
          p_val[i]<-"<0.001"
        }else p_val[i]<-pchisq(wald_weibull[i],1,lower.tail = F)

        }


      AIC=(2*length(fit$par))-2*fit$value
      AICc=AIC-(2*length(fit$par)*(length(fit$par)+1)/(length(t)-length(fit$par)-1))
      BIC=(log(length(t))*length(fit$par))-2*fit$value

      main<-cbind.data.frame(estc,sqrt(diag(varic)),p_val)
      colnames(main)<-c("Estimates","Std Err","P-value")
      rownames(main)<-c("Shape",paste("Shape", co),"Scale",paste("Scale",co),"Intercept",paste("Intercept",co))
      cal<- as.formula(paste("ZIweibull:t ~ ", paste(co, collapse= "+")))

      out<-list(cal,main,fit$value,length(t),sum(t==0),sum(d==0),length(fit$par)-1,AIC,AICc,BIC,wild_boar,fit$hessian)
      names(out)<-c("Call","Coefficients","LogLik","N","Zeros","Censored","df","AIC","AICc","BIC","Data","Hessian")

      return(out)
}else warning('Currently only two variables are possible')

}




} else warning('Only weibull distribution is possible. Make sure dist="weibull" ')


}
