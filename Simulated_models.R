
rm(list=ls(all=TRUE)) 

#Packages to be used
packages<-c("here")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


###################
# Data generation #
###################


source(here("Scripts","Simulation_weibull_generation.R"))
source(here("Scripts","Simulation_gamma_generation.R"))


##################
# ZIWeibull model#
##################


source(here("Scripts","fit_ZIWeibull.R"))



################
# ZIGamma model#
################


source(here("Scripts","fit_ZIGamma.R"))


###################
# Simulation plots#
###################

source(here("Scripts","Graphics_beta.R"))

##################
# Crossed tables #
##################

source(here("Scripts","Crossed_results.R"))