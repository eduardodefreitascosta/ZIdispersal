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

          gamm_zer = (exp(varp[3]))/(1+exp(varp[4]))

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


cbind(estc,sqrt(diag(varic)),estc-1.96*sqrt(diag(varic)),estc+1.96*sqrt(diag(varic)))


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


      cbind(estc,sqrt(diag(varic)),estc-1.96*sqrt(diag(varic)),estc+1.96*sqrt(diag(varic)))
}else warning('Currently only two variables are possible')

}




} else warning('Only weibull distribution is possible. Make sure dist="weibull" ')


}
