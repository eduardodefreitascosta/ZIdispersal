[![DOI](https://zenodo.org/badge/DOI/10.1007/s42081-021-00124-0.svg)](https://link.springer.com/article/10.1007/s42081-021-00124-0)



## This is a paper describying two zero-inflated-right-censored regression models (Weibull and gamma) for wild boar displacement :boar:

## The paper is available online:

https://link.springer.com/article/10.1007/s42081-021-00124-0

de Freitas Costa, E., Schneider, S., Carlotto, G.B. et al. Zero-inflated-censored Weibull and gamma regression models to estimate wild boar population dispersal distance. Jpn J Stat Data Sci (2021). https://doi.org/10.1007/s42081-021-00124-0



# Instrucitons to reproduce the results from the paper "Zero-inflated-censored Weibull and gamma regression models to estimate wild boars population dispersal distance".




### 1 Download the files to a local or cloud folder. The best option is to download all files zipped;

### 2 Unzip and open the R project "ZIdispersal.Rproj" (it is highly recommended to use RStudio);

### 3 In R, you can open the file\ 

  + "Applied_models.R" (available in the tab "Files")\ 
  
  OR\ 
  
  + "Simulated_models.R" (available in the tab "Files")\ The first runs the applied models and the second the simulated.\ 

### 4 Regardless which one (i.e., simulated or applied), first execute the data generation and then the models;

### Note that several folders were created:
  - **Output**: Stores the datasets generated and the results of the simulations models;
  
  - **Figure_applied**: Stores the figures of the simulated ZIWeibull and ZIGamma applied models ;
  
  - **Tables_applied**: Stores the tables generated in the ZIWeibull and ZIGamma applied models;
  
  - **Figure_simulation**: Stores the figures of the simulated ZIGamma model (gamma data) and ZIWeibull model (Weibull data);
  
  - **Tables_simulation**: Stores the tables of the simulated ZIGamma model (gamma data) and ZIWeibull model (Weibull data);
  
  - **Figure_cross**: Stores the figures for both ZIGamma and ZIWeibull using Weibull and gamma data, respectively.
