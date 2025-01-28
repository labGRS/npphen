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
pacman::p_load(tidyverse, terra, sf, ks, tictoc, npphen, scales, microbenchmark)
## Clean environment ----
rm(list = ls(all = T))

## source ExtremeAnoMap -----
source("C:/github/npphen/resolution_test/RFUN_modified/ExtremeAnoMap_nstage1.R")
source('R/ExtremeAnoMap.R')
source('C:/github/npphen/resolution_test/RFUN_modified/ExtremeAnom_nstage1.R')
source('R/Phen.R')
source('C:/github/npphen/resolution_test/RFUN_modified/Phen_extract.R')
## Main not vectorized functions were loaded from npphen main branch
##------------------------------------------------------------------------#
# 3.- General Inputs -----
## Vector benchmarking ----
### loading data ----
f <- system.file("extdata/MegaDrought_spatRast.rda", package = "npphen")
MegaDrought <- readRDS(f)

# Dates
data("modis_dates")

# Generate a Raster time series using a raster stack and a date database from Central Chile
# Obtain data from a particular pixel generating a time series
md_pixel <- cellFromXY(MegaDrought, cbind(313395, 6356610))
md_pixelts <- as.numeric(MegaDrought[md_pixel])

data_test <- tibble(dates = modis_dates, EVI = md_pixelts)

## black and white plot
g <- ggplot(data_test, aes(x = dates, y = EVI/10000, colour = EVI/10000)) +
  geom_line() +
  scale_color_distiller(palette = 'BrBG', 
                        direction = 1) +
  labs(x = '', y = 'EVI', color = 'EVI') +
  scale_x_date(date_breaks = '2 years', date_labels = '%Y') +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.border = element_rect(color = 'white', fill = NA),
        panel.grid.minor = element_blank(),, panel.grid.major = element_blank(), legend.position = 'bottom',
        axis.text = element_text(color = 'white'), axis.title = element_text(color = 'white'), 
        legend.text = element_text(color = 'white'), legend.title = element_text(color = 'white'))
g

ggsave('resolution_test/00_black_white_plot.png', g, width = 10, height = 2.5, dpi = 300)

### timing code ----
#### anomalies ----
times_extreme_anom <- microbenchmark(ExtremeAnom_npphen = val_main <- ExtremeAnom_main(
  x = md_pixelts, dates = modis_dates,
  h = 2, refp = c(1:423), anop = c(1:929), rfd = 0.9, output = "both", rge = c(0, 10000)
),
ExtremeAnom_vectorized = val_vec <- ExtremeAnom(
  x = md_pixelts, dates = modis_dates,
  h = 2, refp = c(1:423), anop = c(1:929), rfd = 0.9, output = "both", rge = c(0, 10000)
),
ExtremeAnom_vec_nstage1 = val_vec_nstage1 <- ExtremeAnom_nstage1(
  x = md_pixelts, dates = modis_dates,
  h = 2, refp = c(1:423), anop = c(1:929), rfd = 0.9, output = "both", rge = c(0, 10000)
), times = 100)

times_plot <- times_extreme_anom %>% autoplot(color = 'grey')

times_plot <- times_plot + 
  theme_minimal() +
  theme(panel.background = element_blank(), panel.border = element_rect(color = 'white', fill = NA),
        panel.grid.minor = element_blank(),, panel.grid.major = element_blank(), legend.position = 'bottom',
        axis.text = element_text(color = 'white'), axis.title = element_text(color = 'white'), 
        legend.text = element_text(color = 'white'), legend.title = element_text(color = 'white')) 

times_plot
ggsave('resolution_test/01_times_plot_ExtremeAnom.png', times_plot, width = 5, height = 3, dpi = 150)


##### Vector differences ----
anom_main <- val_main[str_which(names(val_main), 'anom')]
anom_vec <- val_vec[str_which(names(val_vec), 'anom')]
anom_vec_nstage1 <- val_vec_nstage1[str_which(names(val_vec_nstage1), 'anom')]

rfd_main <- val_main[str_which(names(val_main), 'rfd')]
rfd_vec <- val_vec[str_which(names(val_vec), 'rfd')]
rfd_vec_nstage1 <- val_vec_nstage1[str_which(names(val_vec_nstage1), 'rfd')]

all.equal(anom_main, anom_vec)
all.equal(anom_main, anom_vec_nstage1)

anom_dif <- calcular_diferencias_y_graficar(anom_main, anom_vec, anom_vec_nstage1)
anom_dif

rfd_diff <- calcular_diferencias_y_graficar(rfd_main, rfd_vec, rfd_vec_nstage1)
  
ggsave(anom_dif[[1]], filename = 'resolution_test/02_anom_dif_ExtremeAnom.png', width = 5, height = 3, dpi = 150)
ggsave(rfd_diff[[1]], filename = 'resolution_test/03_rfd_dif_ExtremeAnom.png', width = 5, height = 3, dpi = 150)


#### phenology ----
phen_times <- microbenchmark(
  Phen_npphen = phen_main <- Phen_main(x = md_pixelts, dates = modis_dates, h = 2, frequency = "daily", rge = c(0, 10000)),
  Phen_vectorized = phen_vec <- Phen(x = md_pixelts, dates = modis_dates, h = 2, frequency = "daily", rge = c(0, 10000)),
  Phen_nstage1 = phen_nstage1 <- Phen_nstage1(x = md_pixelts, dates = modis_dates, h = 2, frequency = "daily", rge = c(0, 10000)) , 
  times = 100
)

phen_times_plot <- phen_times %>% autoplot(color = 'grey') +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.border = element_rect(color = 'white', fill = NA),
        panel.grid.minor = element_blank(),, panel.grid.major = element_blank(), legend.position = 'bottom',
        axis.text = element_text(color = 'white'), axis.title = element_text(color = 'white'), 
        legend.text = element_text(color = 'white'), legend.title = element_text(color = 'white'))

phen_times_plot
ggsave('resolution_test/04_times_plot_Phen.png', phen_times_plot, width = 5, height = 3, dpi = 150)

phen_table <- tibble(
  days = 1:length(phen_main),
  main = phen_main,
  vec = phen_vec,
  nstage1 = phen_nstage1
)

phen_plot_series <- phen_table %>% pivot_longer(cols = c('main', 'vec', 'nstage1'), names_to = 'method', values_to = 'value') %>%
  ggplot(aes(x = days, y = value/10000, color = method)) +
  scale_color_manual(values = c('white', 'red', 'cyan'), labels = c('npphen', 'vectorized', 'vectorized nstage1')) +
  labs(x = 'DGS', y = 'EVI', color = '') +
  geom_line()  +
  theme_minimal() +
  theme(panel.background = element_blank(), panel.border = element_rect(color = 'white', fill = NA),
        panel.grid.minor = element_blank(),, panel.grid.major = element_blank(), legend.position = 'bottom',
        axis.text = element_text(color = 'white'), axis.title = element_text(color = 'white'), 
        legend.text = element_text(color = 'white'), legend.title = element_text(color = 'white'))

phen_plot_series
ggsave('resolution_test/05_phen_plot_series.png', phen_plot_series, width = 5, height = 5, dpi = 150)

phen_dif <- calcular_diferencias_y_graficar(phen_main, phen_vec, phen_nstage1)
phen_dif

ggsave(phen_dif[[1]], filename = 'resolution_test/06_phen_dif_Phen.png', width = 5, height = 3, dpi = 150)
##------------------------------------------------------------------------#
# 4.- Raster benchmarking ------
## list files ----
inpath <- 'resolution_test/02_campana_resampled/02a_evi_bands/'

res <- 300
pattern <- paste0('*', res, 'm.tif')
evi_data <- list.files(path = inpath, pattern = glob2rx(pattern), full.names = T)

## dates -----
all_dates <- evi_data %>% str_extract("\\d{4}-\\d{2}-\\d{2}") %>% ymd()

## Anomaly calc ----------
### 300 meters -----
evi_spat <- evi_data %>% rast()
plot(evi_spat[[1:4]])

### Setup ------
nc1 <- 3
refp_end <- which(all_dates <= '2010-06-30') %>% tail(1)
anop_end <- length(all_dates)

outname_npphen <- paste0('resolution_test/02_campana_resampled/trash/anomRFD_campana_2000_2023_res_', res, 'm.tif')
outname_nstage1 <- paste0('resolution_test/02_campana_resampled/trash/anomRFD_campana_2000_2023_res_', res, 'm_nstage1.tif')

tic('ExtremeAnom original flavour')
ExtremeAnoMap(
  s = evi_spat, dates = all_dates, h = 2, refp = c(1:refp_end),
  anop = c(1:anop_end), rfd = 0.9, output = "both", nCluster = nc1, outname = outname_npphen,
  datatype = "INT2S", rge = c(0, 10000)
  )
toc()

tic('ExtremeAnom nstage1')
ExtremeAnoMap_nstage1(
  s = evi_spat, dates = all_dates, h = 2, refp = c(1:refp_end),
  anop = c(1:anop_end), rfd = 0.9, output = "both", nCluster = nc1, outname = outname_nstage1,
  datatype = "INT2S", rge = c(0, 10000)
)
toc()

### 30 meters -----
rm(list = ls(all = T));gc()

inpath <- 'resolution_test/01_campana_30/01a_evi_bands/'
res <- 30

evi_data <- list.files(path = inpath, pattern = glob2rx('*tif'), full.names = T)

### dates -----
all_dates <- evi_data %>% str_extract("\\d{4}-\\d{2}-\\d{2}") %>% ymd()



evi_spat <- evi_data %>% rast()
plot(evi_spat[[1:4]])

### Setup ------
nc1 <- 2
refp_end <- which(all_dates <= '2010-06-30') %>% tail(1)
anop_end <- length(all_dates)

outname_npphen <- paste0('resolution_test/01_campana_30/trash/anomRFD_campana_2000_2023_res_', res, 'm.tif')
outname_nstage1 <- paste0('resolution_test/01_campana_30/trash/anomRFD_campana_2000_2023_res_', res, 'm_nstage1.tif')

tic('ExtremeAnom original flavour 30')
ExtremeAnoMap(
  s = evi_spat, dates = all_dates, h = 2, refp = c(1:refp_end),
  anop = c(1:anop_end), rfd = 0.9, output = "both", nCluster = nc1, outname = outname_npphen,
  datatype = "INT2S", rge = c(0, 10000)
)
toc()

tic('ExtremeAnom nstage1 30')
ExtremeAnoMap_nstage1(
  s = evi_spat, dates = all_dates, h = 2, refp = c(1:refp_end),
  anop = c(1:anop_end), rfd = 0.9, output = "both", nCluster = nc1, outname = outname_nstage1,
  datatype = "INT2S", rge = c(0, 10000)
)
toc()

##------------------------------------------------------------------------#
# 5.- check raster differences ----
rm(list = ls(all = T)); gc()

## 300 meters ----
inpath_anom <- 'resolution_test/02_campana_resampled/02b_anom_bands/'
inpath_rfd <- 'resolution_test/02_campana_resampled/02c_rfd_bands/'

### list files
anom_vec_300_ls <- list.files(path = inpath_anom, pattern = glob2rx('*300m.tif'), full.names = T)
anom_nstage_300_ls <- list.files(path = inpath_anom, pattern = glob2rx('*nstage1.tif'), full.names = T)

rfd_vec_300_ls <- list.files(path = inpath_rfd, pattern = glob2rx('*300m.tif'), full.names = T)
rfd_nstage_300_ls <- list.files(path = inpath_rfd, pattern = glob2rx('*nstage1.tif'), full.names = T)

## read anom ----
anom_rfd_300_vec <- rast(anom_vec_300_ls)
anom_rfd_300_nstage1 <- rast(anom_nstage_300_ls)

### diff anom ------
difference <- function(x, y){
  if(all(is.na(x)) | all(is.na(y))){
    return(rep(NA, length(x)))
  }

  diff <- abs(y) - abs(x)
  return(diff)
}

diff_anom_300 <- app(c(anom_rfd_300_vec, anom_rfd_300_nstage1), fun = function(a){
  idx <- length(anom_vec_300_ls)+1
  idx_end <- length(anom_vec_300_ls)*2
  x <- a[1:length(anom_nstage_300_ls)]
  y <- a[idx:idx_end]
  
  difference(x, y)
})



plot(diff_anom_300[[1:4]])

mean_diff <- diff_anom_300 %>% mean(na.rm = T)
sd_diff <- diff_anom_300 %>% app(fun = sd, na.rm = T) 

writeRaster(c(mean_diff, sd_diff), 'resolution_test/02_campana_resampled/02e_diff_comparison/anomFunction_diff_compare_300.tif', datatype = 'INT2S', overwrite = T)

## read rfd -----
rfd_300_vec <- rast(rfd_vec_300_ls)
rfd_300_nstage1 <- rast(rfd_nstage_300_ls)

### diff rfd ------  

diff_rfd_300 <- app(c(rfd_300_vec, rfd_300_nstage1), fun = function(a){
  idx <- length(anom_vec_300_ls)+1
  idx_end <- length(anom_vec_300_ls)*2
  x <- a[1:length(anom_nstage_300_ls)]
  y <- a[idx:idx_end]
  
  difference(x, y)
})

plot(diff_rfd_300[[1:4]])

mean_diff <- diff_rfd_300 %>% mean(na.rm = T)
sd_diff <- diff_rfd_300 %>% app(fun = sd, na.rm = T) 

writeRaster(c(mean_diff, sd_diff), 'resolution_test/02_campana_resampled/02e_diff_comparison/rfdFunction_diff_compare_300.tif', datatype = 'INT2S', overwrite = T)
