
model <- function(N,varp){
  a       =  exp(varp[1])*exp(y*varp[2])
  lam     =  exp(varp[3])*exp(y*varp[4])
  
  # Parameters - WEIGHT IN ZEROS AND INFS: GAMMA ZERO AND GAMMA ONES WITH A COVARIABLE
 # gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6]))  
 # gamm_inf = (exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  
  U=runif(N)
  V=runif(N,gamm_zer,1)
  y=numeric(N)
  for(i in 1:N){
    if(U[i]<=gamm_zer[i]){ 
      y[i] = 0  
    }else 
      if(U[i]>gamm_zer[i] & U[i] <= 1){
        #y[i]=qweibull(V[i], a, lam)
        y[i] = lam[i]*((-log(1 - (V[i]-gamm_zer[i])/(1-gamm_zer[i])))^(1/a[i]))
      }else {      
        y[i]=+Inf
      }
  }
  return(y)
}


GPE=function(varp)
{
  a       =   exp(varp[1])*exp(y*varp[2])
  lam     =   exp(varp[3])*exp(y*varp[4])
  #gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
  gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6]))
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