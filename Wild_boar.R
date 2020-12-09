
wild_boar<-function(){
#Sample size
n <- 400

set.seed(13)
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
wild_boar <- as.data.frame(cbind(t,Z,delta, ages,x1))
names(wild_boar) <- c("dist","zero","delta","age", "sex")

saveRDS(wild_boar,file = here("wild_boar.rds"))
}
