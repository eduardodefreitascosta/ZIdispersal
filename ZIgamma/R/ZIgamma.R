

ZIgamma<-function(ze,t,d,data,co=c(),dist){

  #Packages to be used
  packages<-c("flexsurv","survival","here")


  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))


  attach(data)

  if(dist=="gamma"){

    if(length(co)!=0){

      data_log<-cbind.data.frame(ze,data[, c(co)])


      a<-glm(ze~.,family=binomial(link='logit'),data=data_log)

      (fmla <- as.formula(paste("Surv(t,d,type='right') ~ ", paste(co, collapse= "+"))))

      data_gamma<-(cbind.data.frame(t,d, data_log))

      b<-flexsurvreg(fmla,data=data_gamma,subset=t!=0,dist=dist)

    }else{


      a<-glm(ze~1,family=binomial(link='logit'))

      (fmla <- paste("Surv(t,d,type='right') ~ 1"))


      b<-flexsurvreg(Surv(subset(t,t>0),subset(d,t>0),type='right') ~ 1,dist=dist)

    }


    summa<-list(a,b)
    names(summa)<-c("Logistic","Gamma")

    detach(data,unload = T)
    return(summa)
  }else warning('So far, only gamma distribution is possible:"gamma" ')

}



ZIgamma(ze=zero,t=dist,d=delta,data=wild_boar,co=c("age"),dist="weibull")
