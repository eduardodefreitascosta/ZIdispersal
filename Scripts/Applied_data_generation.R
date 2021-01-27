
#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","flexsurv","knitr","glmmsr","plotly","gridExtra","grid","ggridges","ggthemes","ggcorrplot","survival")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))



wd <- getwd()
set.seed(13)



#Sample size
n <- 400


#Creating sex
z<- rep(1,n/2) 
w<-rep(0,n/2)
x1<-  c(z,w)

#Creating age in monts
age_fem<-rexp(n/2,1/12)
age_mac<-rexp(n/2,1/15)
ages<-c(age_fem,age_mac)

# Crianting distances, zero and censure. 
zer<- (exp(-1.5 - 0.02*(ages) + 0.2*x1))/(1+exp(-1.5 - 0.02*(ages) + 0.2*x1))  #Probabilidade de ser zero


alpha.t <- 1 #makes exponential
lambda.t <- 0.1

#fem
gamma.t <- lambda.t*exp(-0.02*(age_fem)+ 0.5*z)
u <- runif(n)
y_fem_auxi <- NULL

for(i in 1:(n/2)){
  y_fem_auxi[i] <- (-log(u[i])/gamma.t[i])^(1/alpha.t) 
}

y_fem <- y_fem_auxi*(1-rbinom(200,1,zer[1:200]))  #Distances for females

summary(y_fem)
hist(y_fem)

#mac
gamma.t <- lambda.t*exp(-0.02*(age_mac) + 0.5*w)
u <- runif(n)
y_mac_auxi <- NULL

for(i in 1:(n/2)){
  y_mac_auxi[i] <- (-log(u[i])/gamma.t[i])^(1/alpha.t) 
}

y_mac <- y_mac_auxi*(1-rbinom(200,1,zer[201:400])) #Distances for males

summary(y_mac)
hist(y_mac)

y<-c(y_fem, y_mac ) 

c <- rep(30,n)

t <- NULL
delta <- NULL
Z <- NULL
for (i in 1:n){ 
  t[i] = min(y[i],c[i]) #observed distance (minimum between failure and cens)
  
  if (t[i] < c[i]){
    delta[i]=1          #delta: indicative of failures and cens
  } else {delta[i]=0}
  
  if (t[i] == 0){ #Z: zero inticative
    Z[i]=1
  }else{Z[i]=0}
}

#Creating dataframe
data <- as.data.frame(cbind(t,Z,delta, ages,x1))
names(data) <- c("dist","zero","delta","age", "sex")  

#percentage of events
sum(data$delta)/n
#[1] 0.9725

#percentage of zeros
sum(data$zero)/n
#[1] 0.2375

plot(Surv(data$dist,delta))


##################
#Data summary    #
##################

dir.create(paste(wd,"/Tables_applied",sep=""))

tabela<-cbind(
  #General mean distande and per sex
  rbind(
    mean(data$dist[data$sex==0]), # geral para machos
    mean(data$dist[data$sex==0&data$age<15]), #Machos <16
    mean(data$dist[data$sex==0&data$age>15]), #Machos >16
    mean(data$dist[data$sex==1]), # geral para femeas
    mean(data$dist[data$sex==1&data$age<12]), #Femeas <16
    mean(data$dist[data$sex==1&data$age>12]) #Femeas >16
  ),
  
  # %zeros general mean and per sex
  rbind(
    mean(data$zero[data$sex==0]), # general for males
    mean(data$zero[data$sex==0&data$age<15]), #males <16
    mean(data$zero[data$sex==0&data$age>15]), #males >16
    mean(data$zero[data$sex==1]), # general for females
    mean(data$zero[data$sex==1&data$age<12]), #females <16
    mean(data$zero[data$sex==1&data$age>12]) #females >16
  ),
  
  
  #% General mean censoriing and per sex
  rbind(
    mean(data$delta[data$sex==0]), # general for males
    mean(data$delta[data$sex==0&data$age<15]), #males <16
    mean(data$delta[data$sex==0&data$age>15]), #males >16
    mean(data$delta[data$sex==1]), # general for females
    mean(data$delta[data$sex==1&data$age<12]), #Females <16
    mean(data$delta[data$sex==1&data$age>12]) #Females >16
  )
  
)

row.names(tabela)<-c("Males","Males <16","Males >16","Females","Females <16","Females >16")
colnames(tabela)<-c("Mean dispersal","%zeros","%Recapture")
kable(tabela)

write.table(tabela,here("Tables_applied","general.txt"),sep=";")


dir.create(paste(wd,"/Figure_applied",sep=""))

data$Sex<-c(rep("Female",200),rep("Male",200))

jpeg(here("Figure_applied","dist_hist.jpg"), height = 4, width = 6, units = 'in', res=300)
ggplot(data, aes(x=dist))+
  geom_histogram(color="black", fill="white")+
  facet_grid(Sex ~ .)+
  labs(x ="Distance (Km)", y = "Frequency")
dev.off()
