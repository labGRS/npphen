# 1.- Phenology experiments: Extraction and C++ optimization ----------
##------------------------------------------------------------------------#
## @date 2025-01-20
## @project C:/github/npphen
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
##------------------------------------------------------------------------#
# 2.- Libraries --------
pacman:p_load(tidyverse, terra, sf, ks)
## Clean environment ----
rm(list = ls(all = T))

##------------------------------------------------------------------------#
# 3.- General Inputs -----
## base Phen function ------
source('R/Phen.R')

## Sample time series -----
data("phents")
all_dates <- as.Date(phents$dates)

