

ZIgamma<-function(ze,t,d,data,co=c()){


  data_log<-data[, c(co)]

  #Packages to be used
  packages<-c("flexsurv","survival","here")


  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  wild_boar <- readRDS(here("ZIgamma","wild_boar.rds"))

  a<-glm(ze~.,family=binomial(link='logit'),data=data_log)

  (fmla <- as.formula(paste("Surv(t,d,type='right') ~ ", paste(co, collapse= "+"))))

  data_gamma<-(cbind.data.frame(t,d, data_log))

  b<-flexsurvreg(fmla,data=data_gamma,subset=t!=0,dist="gamma")



  summa<-list(a,b)
  names(summa)<-c("Logistic","Gamma")


  return(summa)


}



ZIgamma(ze=zero,t=dist,d=delta,data=wild_boar,co=c("age","sex"))
