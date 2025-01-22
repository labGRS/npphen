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
pacman::p_load(tidyverse, terra, sf, ks, parallel)
## Clean environment ----
rm(list = ls(all = T))

##------------------------------------------------------------------------#
# 3.- General Inputs -----
## base Phen function ------
ls_f <- list.files(path = "R/", pattern = "*.R", full.names = T)
sapply(ls_f, source)

## Sample time series -----
data("phents")
all_dates <- as.Date(phents$dates)

## Benchmarking ----
times <- microbenchmark::microbenchmark(
  time_phen_vectorized = Phen(x = phents$NDVI, dates = all_dates, h = 2, frequency = "daily", rge = c(0, 10000), plot = F),
  time_phen_original = npphen::Phen(x = phents$NDVI, dates = all_dates, h = 2, frequency = "daily", rge = c(0, 10000), plot = F)
)
# Median 3.98   # min 3.4    # lq 3.72
# mean 4.37     # max 7.74   # uq 4.67
autoplot(times)

## Benchmarking ----
### HPI ----
hpi_time <- microbenchmark::microbenchmark(
  hpi = Hmat <- ks::Hpi(na.omit(D1)), ## using nstage = 1 seems to improve speed without too much compromise of precision
  hpi_nstage_1 = Hmat_2 <- ks::Hpi(na.omit(D1), nstage = 1)
)
autoplot(hpi_time)
## Median ~2 seconds

### KDE ----
kde_time <- microbenchmark::microbenchmark(
  kde = K1 <- ks::kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500)),
  kde_time = K1_2 <- ks::kde(na.omit(D1), H = Hmat_2, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
)

autoplot(kde_time)
# median ~ 3.8

##------------------------------------------------------------------------#
# 4.- Extract phen  ----
## source modified -----
source("resolution_test/RFUN_modified//Phen_extract.R")
## initial arguments ----
x <- phents$NDVI
dates <- all_dates
h <- 2
frequency <- "daily"
rge <- c(0, 10000)
plot <- F
anop <- 1:length(dates)
refp <- 1:length(dates)
output <- 'both'

## test_data -----
f <- system.file("extdata/Bdesert_spatRast.rda", package = "npphen")
test <- readRDS(f)

### 1 pixel ----
px <- cellFromXY(test, cbind(286638, 6852107))
px_ts <- as.numeric(test[px])

times_phen <- microbenchmark::microbenchmark(
  phen_02 =  px_val_02 <- Phen_extract(x = px_ts, dates = dates, h = h, frequency = frequency, rge = rge, plot = plot),
  phen_01 =  px_val_01 <- Phen(x = px_ts, dates = dates, h = h, frequency = frequency, rge = rge, plot = plot),
  times = 10
) 
autoplot(times_phen)

## compare values
plot(px_val_01, type = 'l')
lines(px_val_02[[1]], col = 'red')


### map example with CumDensity -----
cores <- 6
cluster <- makeCluster(cores)
clusterExport(cl = cluster, varlist = c("Phen_extract", "dates", "h", "frequency", "rge", "plot"))

map_out <- app(test, fun = function(x) {
  Phen_extract(x = x, dates = dates, h = h, frequency = frequency, rge = rge, plot = plot)
}, cores = cluster)

stopCluster(cluster)
