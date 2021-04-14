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




###--------------------------------------------------------------------------------------------###
# Graphics for gamma distribution
###--------------------------------------------------------------------------------------------###

# Define files
#dir.create(paste(wd,"/Figure_gamma",sep=""))


# Set environment
lime = 10
N<-c(200,400,600,800,1000)
modelos <- c("dados_No_p1", "modelos_No_p1")
samp<-c(200,400,600,800,1000)

##Remember the parameters
beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1
varp=c(beta1,beta2,beta3,beta4,beta5,beta6)

#####Betas ZIWeibull model
betas200 <- read.table(paste(wd,'/','/Output/simulation_gamma',N[1],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas400 <- read.table(paste(wd,'/','/Output/simulation_gamma',N[2],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas600 <- read.table(paste(wd,'/','/Output/simulation_gamma',N[3],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas800 <- read.table(paste(wd,'/','/Output/simulation_gamma',N[4],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))
betas1000 <- read.table(paste(wd,'/','/Output/simulation_gamma',N[5],'/','modelos_No_p1/weibull/param_Weibull.txt',sep=''))

betas<-cbind(apply(betas200[,1:6],2,mean),
             apply(betas400[,1:6],2,mean),
             apply(betas600[,1:6],2,mean),
             apply(betas800[,1:6],2,mean),
             apply(betas1000[,1:6],2,mean))
vies<-(betas-varp)

###########Betas ZIGamma model
gammas200 <- read.table(paste(wd,'/','Output/simulation',N[1],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas400 <- read.table(paste(wd,'/','Output/simulation',N[2],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas600 <- read.table(paste(wd,'/','Output/simulation',N[3],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas800 <- read.table(paste(wd,'/','Output/simulation',N[4],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))
gammas1000 <- read.table(paste(wd,'/','Output/simulation',N[5],'/','modelos_No_p1/gamma/param_gamma.txt',sep=''))



#Transform the parameters to natural scale
#V1: beta5=-3
#V2: beta6=1
#V3: beta1=0.5
#V4: beta3=1.5
#V5: beta4=2
#V6: beta2=0.5
#V7: LogLik
#V7: AIC


gammas<-cbind(apply(gammas200,2,mean),
              apply(gammas400,2,mean),
              apply(gammas600,2,mean),
              apply(gammas800,2,mean),
              apply(gammas1000,2,mean)
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
dps_g<-cbind(apply(gammas200,2,sd),
             apply(gammas400,2,sd),
             apply(gammas600,2,sd),
             apply(gammas800,2,sd),
             apply(gammas1000,2,sd))


variancia_g<-cbind(apply(gammas200,2,var),
                   apply(gammas400,2,var),
                   apply(gammas600,2,var),
                   apply(gammas800,2,var),
                   apply(gammas1000,2,var))

# Standard errors
#Weibull model
var200 <- read.table(paste(wd,"/Output/simulation_gamma200/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var400 <- read.table(paste(wd,"/Output/simulation_gamma400/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var600 <- read.table(paste(wd,"/Output/simulation_gamma600/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var800 <- read.table(paste(wd,"/Output/simulation_gamma800/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))
var1000 <- read.table(paste(wd,"/Output/simulation_gamma1000/modelos_No_p1/weibull/Erro_Weibull.txt",sep=""))



vars<-cbind(apply(var200[,1:6],2,mean),
            apply(var400[,1:6],2,mean),
            apply(var600[,1:6],2,mean),
            apply(var800[,1:6],2,mean),
            apply(var1000[,1:6],2,mean))

eqm<-(vies^2+variancia)

#Gamma model
varg200 <- read.table(paste(wd,"/Output/simulation200/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg400 <- read.table(paste(wd,"/Output/simulation400/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg600 <- read.table(paste(wd,"/Output/simulation600/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg800 <- read.table(paste(wd,"/Output/simulation800/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))
varg1000 <- read.table(paste(wd,"/Output/simulation1000/modelos_No_p1/gamma/Erro_gamma.txt",sep=""))

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




###--------------------------------------------------------------------------------------------###
# Tables for crossed data and model
###--------------------------------------------------------------------------------------------###

# Define files

dir.create(paste(wd,"/Tables_cross",sep=""))


## Theoretic E(t|x=) for gamma
media1<-function(x){
  alfa=exp(beta1+beta2*x)
  teta=exp(beta3+beta4*x)
  p0=exp(beta5+beta6*x)/(1+exp(beta5+beta6*x))
  
  media.1=(1-p0)*teta*alfa
  return(rbind(media.1,p0))
}

media1(-1);media1(1)




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

means_w

colnames(means_w)[1:5]<-c("200","400","600","800","1000")

row.names(means_w)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(means_w,here("Tables_cross","means_weibull.txt"),sep=";",row.names = T)


## Theoretic E(t|x=) for Weibull
media<-function(x){
  alfa=exp(beta1+beta2*x)
  teta=exp(beta3+beta4*x)
  p0=exp(beta5+beta6*x)/(1+exp(beta5+beta6*x))
  
  media.1=(1-p0)*teta*gamma(1+(1/alfa))
  return(rbind(media.1,p0))
}

media(-1);media(1)


##Mean time ZIGamma 
po<-NULL
scale<-NULL
shape<-NULL

mediag<-function(x){
  for (j in 1:5){
    po[j]=exp(gammas[5,j]+gammas[6,j]*x)/(1+exp(gammas[5,j]+gammas[6,j]*x))
    scale[j]=exp(gammas[3,j]+gammas[4,j]*x)
    shape[j]=exp(gammas[1,j]+gammas[2,j]*x)
  }
  media.g=(1-po)*(shape)*(scale)
  return(rbind(media.g,po))
}

means_g<-rbind(
  mediag(1),
  mediag(-1)
)

means_g

colnames(means_g)[1:5]<-c("200","400","600","800","1000")

row.names(means_g)<-c("Mean x=1","P0 x=1","Mean x=-1","P0 x=-1")

write.table(round(means_g,3),here("Tables_cross","means_gamma.txt"),sep=";",row.names = T)


