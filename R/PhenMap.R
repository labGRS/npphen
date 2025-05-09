#' @title PhenMap
#' @description Estimates annual Land Surface Phenology (LSP) using time series of a vegetation greenness raster stack.
#' @encoding UTF-8
#' @param s SpatRaster of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior. The code has been optimized to work with integer values. Please re-scale the input SpatRaster if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000)
#' @param dates A date vector. The number of dates must be equal to the number of layers of the vegetation index SpatRaster.
#' @param h Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h = 1 if the vegetation is in the Northern Hemisphere (season starting on January 1st), h = 2 if it is in the Southern Hemisphere (season starting on July 1st)
#' @param frequency Character string. Defines the number of samples for the output phenology and must be one of the following options: 'daily' with an output vector of length 365, '8-days' with an output vector of length 46 (e.g. MOD13Q1 and MYD13Q1 combined), 'monthly' with an output vector of length 12,'bi-weekly' with an output vector of length 24 (i.e. GIMMS) or '16-days' (default) with an output vector of length 23 (i.e MOD13Q1 or MYD13Q1)
#' @param nCluster Numeric. Number of CPU cores to be used for computational calculations
#' @param outname Character vector with the output directory path and filename with extension or only the filename and extension if working directory was set. For example outname = "output_phen.tif". See \code{\link[terra]{writeRaster}}
#' @param datatype Character. Output data type. For example datatype = "INT2S" for signed integer values. See \code{\link[terra]{writeRaster}}
#' @param rge  A vector containing minimum and maximum values of the vegetation index used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000)
#' @details Per pixel, it derives the most recurrent annual Land Surface Phenology (LSP) for a reference period  using time series of satellite based greenness values such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). The annual LSP is calculated for all pixels of the input SpatRaster in the same way as for the \code{\link{Phen}} function. The output is a multiband raster where every band is the expected greenness value at a given time step of the growing season. For example, for MODIS Vegetation Index 16-days composites the number of time steps of the growing season is 23 (frequency = "16-days"), and therefore, the output raster will have 23 bands. A vector with dates for the greenness values is also required
#' @return SpatRaster
#' @seealso \code{\link{Phen}}
#' @examples
#' \dontrun{
#' ## DEPENDING ON HARDWARE, THIS PROCESS CAN BE HIGHLY TIME CONSUMING##
#'
#' ## Testing with an NDVI spatRast from Central Chile, h = 2##
#'
#' # Load data
#' # SpatRaster
#' library(terra)
#' f <- system.file("extdata/MegaDrought_spatRast.rda", package = "npphen")
#' MegaDrought <- readRDS(f)
#' # Dates
#' data("modis_dates")
#'
#' # Making the LSP output raster, n bands = 23
#'
#' # Define the number of cores to be use. In this example we use 1
#' nc1 <- 1
#'
#' PhenMap(
#'   s = MegaDrought, dates = modis_dates, h = 2,
#'   frequency = "16-days", nCluster = nc1, outname = "phen_MD.tif",
#'   datatype = "INT2S", rge = c(0, 10000)
#' )
#' # map1 <- rast("phen_MD.tif")#run only for load phenology brick
#' # plot(map1)
#'
#' ## Testing with an NDVI spatRast from the Atacama Desert, Northern Chile, h=2 ##
#'
#' # Load data
#' SpatRaster
#' f <- system.file("extdata/Bdesert_spatRast.rda", package = "npphen")
#' Bdesert <- readRDS(f)
#' # Dates
#' data("modis_dates")
#'
#' # Making the LSP output raster, n bands = 23
#' # Define the number of cores to be use. In this example we use 1
#' nc1 <- 1
#'
#' PhenMap(
#'   s = Bdesert, dates = modis_dates, h = 2,
#'   frequency = "16-days", nCluster = 1, outname = "phen_BD.tif",
#'   datatype = "INT2S", rge = c(0, 10000)
#' )
#' # map2 <- rast("phen_BD.tif") #run only for loading the multiband LSP spatRaster
#' # plot(map2)
#' }
#' @export

PhenMap <-
  function(s, dates, h, frequency = "16-days", nCluster, outname, datatype, rge) {
    ff <- function(x) {
      if (length(rge) != 2) {
        stop("rge must be a vector of length 2")
      }
      if (rge[1] > rge[2]) {
        stop("rge vector order must be minimum/maximum")
      }
      if (length(dates) != length(x)) {
        stop("N of dates and files do not match")
      }

      freq.method <- match(frequency, c(
        "daily", "8-days", "16-days",
        "bi-weekly", "monthly"
      ))
      if (is.na(freq.method)) {
        stop("Invalid frequency. Must be one of: daily, 8-days, 16-days, bi-weekly or monthly")
      }
      nGS <- switch(frequency,
        "daily" = 365L,
        "8-days" = 46L,
        "16-days" = 23L,
        "monthly" = 12L,
        "bi-monthly" = 24L
      )
      if (all(is.na(x))) {
        return(rep(NA, nGS))
      }
      if (all(x < rge[1], na.rm = T) || all(x > rge[2], na.rm = T)) {
        return(rep(NA, nGS))
      }

      DOY <- lubridate::yday(dates)
      DOY[which(DOY == 366)] <- 365L
      D1 <- cbind(DOY, x)
      if (length(unique(D1[, 2])) < 10 || (nrow(D1) - sum(is.na(D1))) < (0.1 * nrow(D1))) {
        return(rep(NA, nGS))
      }

      if (h != 1 && h != 2) {
        stop("Invalid h")
      }
      DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
      if (h == 2) {
          D1[, 1] <- DOGS[match(D1[, 1], DOGS[, 1]), 2]
      }

      Hmat <- ks::Hpi(na.omit(D1))
      Hmat[1, 2] <- Hmat[2, 1]
      K1 <- ks::kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
      K1Con <- K1$estimate
      K1Con <- apply(K1$estimate, 1, function(row) {
        Kdiv <- sum(row)
        if (Kdiv == 0) {
          return(rep(0, length(row)))
        } else {
          return(row / Kdiv)
        }
      })

      K1Con <- t(K1Con)

      first.no.NA.DOY <- min(D1[, 1][which(is.na(D1[, 2]) == FALSE)])
      last.no.NA.DOY <- max(D1[, 1][which(is.na(D1[, 2]) == FALSE)])

      n.selects <- apply(K1Con, 1, function(row) which(row == max(row)))

      MAXY <- mapply(function(n.select, i) {
        if (length(n.select) > 1) {
          return(NA)
        } else if (length(n.select) == 1 && i >= first.no.NA.DOY && i <= last.no.NA.DOY) {
          return(median(K1$eval.points[[2]][n.select]))
        } else {
          return(NA)
        }
      }, n.select = n.selects, i = 1L:365L)

      select_DGS <- switch(frequency,
        "daily" = 1L:365,
        "8-days" = seq(1L, 365L, 8L),
        "16-days" = seq(1L, 365L, 16L),
        "monthly" = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
        "bi-weekly" = c(1, 15, 32, 46, 60, 74, 91, 105, 121, 135, 152, 166, 182, 196, 213, 227, 244, 258, 274, 288, 305, 319, 335, 349)
      )

      Ref <- MAXY[select_DGS]

      names(Ref) <- select_DGS
      Ref
    }
    app(s, fun = ff, filename = outname, cores = nCluster, overwrite = T, wopt = list(datatype = datatype))
  }
