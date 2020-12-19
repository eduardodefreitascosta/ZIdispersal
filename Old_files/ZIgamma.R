
ZIgamma<-function(formula,zero,censor,data,dist){

  #Packages to be used
  packages<-c("flexsurv","survival")


  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))

  attach(data)


  data1<-model.frame(formula,data=data)
  mt <- attr(x = data1, which = "terms")

  x <- model.matrix(object = mt, data = data1,xlev = T)

  t<-data1[,1]

  (fmla <- as.formula(paste("zero ~ " , (Reduce(paste,deparse(mt[[3]]))))))

  a<-glm(fmla,family=binomial(link='logit'),data=data)


  (fmla2 <- as.formula(paste(paste("Surv(t,censor,type='right') ~ " , (Reduce(paste,deparse(mt[[3]])))))))


  b<-flexsurvreg(fmla2,data=data,subset=t!=0,dist=dist)


  summa<-list(a,b)
  names(summa)<-c("Logistic","Gamma")

  detach(data,unload = T)
  return(summa)

  }


ZIgamma(dist~sex+age,zero=zero,censor=delta,data=wild_boar,dist="gamma")








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

