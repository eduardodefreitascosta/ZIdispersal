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



################
#Applied models#
################



##Data generation

source(here("Scripts","Applied_data_generation.R"))



##Model application

source(here("Scripts","fit_Applied.R"))
