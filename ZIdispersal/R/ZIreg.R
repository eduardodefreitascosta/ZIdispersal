

ZIreg<-function(formula,zero,censor,data,dist){

attach(data)
  data1<-model.frame(formula,data=data)
  mt <- attr(x = data1, which = "terms")

  x <- model.matrix(object = mt, data = data1,xlev = T)

  t<-data1[,1]

  if(dist=="weibull"){

  l1<-1+length(x[1,]);l2<-length(x[1,])*2;l3<-1+length(x[1,])*2;l4<-length(x[1,])*3


      GPE=function(varp){
          var1<-varp[1:length(x[1,])]
          a       =   exp(x%*%var1)

          var2<-varp[l1:l2]
          lam     =   exp(x%*%var2)

          var3<-varp[l3:l4]
          gamm_zer = exp(x%*%var3) /(1+exp(x%*%var3))

          f0   = gamm_zer
          fw = (exp(log(a)-log(lam)))*((exp(log(t)-log(lam)))^(a-1))*(exp(-(exp(log(t)-log(lam)))^a))
          Sw = (exp(-(exp(log(t)-log(lam)))^a))


        fpop = (1-gamm_zer)*fw
        Spop = (1-gamm_zer)*Sw
        f    = (fpop^censor)*(Spop^(1-censor))
        g   <- ifelse(t==0, f0, f)
        adFunc   = sum(log(g))
        return(adFunc)
      }


      # initial value
      varp=rep(0,l4)

      fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))

      estc = try(fit$par)
      Hc=try(fit$hessian); varic=try(-solve(Hc))

      detach(data,unload = T)


      #Hypothesis test weibull
      betas_weibull=estc
      SE_weibull=sqrt(diag(varic))
      wald_weibull=betas_weibull^2/SE_weibull^2


      p_val<-numeric()
      for (i in 1:length(fit$par)){

        p_val[i]<-round(pchisq(wald_weibull[i],1,lower.tail = F),3)

        if(p_val[i]<0.001){
          p_val[i]<-"<0.001"
        }else p_val[i]<-p_val[i]
      }


      AIC=(2*length(fit$par))-2*fit$value
      AICc=AIC-(2*length(fit$par)*(length(fit$par)+1)/(length(t)-length(fit$par)-1))
      BIC=(log(length(t))*length(fit$par))-2*fit$value

      main<-cbind.data.frame(round(estc,3),round(sqrt(diag(varic)),3),p_val)
      colnames(main)<-c("Estimates","Std Err","P-value")
      rownames(main)<-make.names(c("Shape",names(x[1,])[-1],"Scale",names(x[1,])[-1],"P0",names(x[1,])[-1]), unique=TRUE)



      #cal<- (paste(formula))


      out<-list(print(formula),main,fit$value,length(t),sum(t==0),sum(censor==0),length(fit$par)-1,AIC,AICc,BIC)
      names(out)<-c("Call","Coefficients","LogLik","N","Zeros","Censored","df","AIC","AICc","BIC")

      return(out)

}
  if(dist=="gamma"){

  #Packages to be used
  packages<-c("flexsurv","survival")


  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))


  (fmla <- as.formula(paste("zero ~ " , (Reduce(paste,deparse(mt[[3]]))))))

  a<-glm(fmla,family=binomial(link='logit'),data=data)


  (fmla2 <- as.formula(paste(paste("Surv(t,censor,type='right') ~ " , (Reduce(paste,deparse(mt[[3]])))))))
  
  (fmla3 <- as.formula(paste(paste(" ~ " , (Reduce(paste,deparse(mt[[3]])))))))
  
  custom.gamma <- list(name = dist,
                       pars = c("shape", "scale"),
                       location = "scale",
                       transforms = c(log, log),
                       inv.transforms = c(exp, exp),
                       inits = function(t){
                         c(0.01, 100)},
                       method="Nelder-Mead"
  )
  
  b<-flexsurvreg(fmla2, anc = list(shape = fmla3),
                data=data,subset=t!=0,dist=custom.gamma)


  summa<-list(a,b)
  names(summa)<-c("Logistic","Gamma")

  detach(data,unload = T)
  return(summa)

}else warning('Only weibull and gamma distributions are possible. Make sure dist="weibull" or dist="gamma" ')



}
