install.packages("devtools")
library(devtools)


install_github("eduardodefreitascosta/Zidispersal/ZIdispersal",force=T)
library(ZIdispersal)

ZIreg(dist~age:sex,zero=zero,censor=delta,data=wild_boar,dist="weibull")


