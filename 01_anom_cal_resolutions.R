# 1.- Anomaly calc at different resolutions ----------
##------------------------------------------------------------------------#
## @date 2025-01-20
## @project C:/github/npphen
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
##------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(tidyverse, terra, sf, ks)
## Clean environment ----
rm(list = ls(all = T))

##------------------------------------------------------------------------#
# 3.- General Inputs -----
## source ExtremeAnoMap -----
source("C:/github/npphen/R/ExtremeAnoMap.R")

## list files ----
inpath <- 'resolution_test/02_campana_resampled/02a_evi_bands/'

res <- 0.0009
pattern <- paste0('*res_', res, 'm*.tif')
evi_data <- list.files(path = inpath, pattern = glob2rx(pattern), full.names = T)

## dates -----
all_dates <- evi_data %>% str_extract("\\d{4}-\\d{2}-\\d{2}") %>% ymd()

# 4.- Anomaly calc ----------
evi_spat <- evi_data %>% rast()
plot(evi_spat[[1]])

## Setup ------
nc1 <- 3
refp_end <- which(all_dates <= '2010-06-30') %>% tail(1)
anop_end <- length(all_dates)

outname <- paste0('resolution_test/02_campana_resampled/trash/anomRFD_campana_2000_2023_res_', res, 'm.tif')

ExtremeAnoMap(
  s = evi_spat, dates = all_dates, h = 2, refp = c(1:refp_end),
  anop = c(1:anop_end), rfd = 0.9, output = "both", nCluster = nc1, outname = outname,
  datatype = "INT2S", rge = c(0, 10000)
  )
