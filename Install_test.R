install.packages("devtools")
library(devtools)


install_github("eduardodefreitascosta/Zidispersal/ZIdispersal",force=T)
library(ZIdispersal)

ZIreg(ze=zero,t=dist,d=delta,data=wild_boar,co=c(),dist="gamma")

