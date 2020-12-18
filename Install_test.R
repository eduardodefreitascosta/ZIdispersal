install.packages("devtools")
library(devtools)


install_github("eduardodefreitascosta/Zidispersal/ZIdispersal",force=T)
library(ZIdispersal)

ZIreg(distãge+sex,zero=zero,censor=delta,data=wild_boar,dist="gamma")

