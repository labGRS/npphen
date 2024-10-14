#' @title PhenKplot
#' @description Plot the most probable vegetation greenness values.
#' @encoding UTF-8
#' @param x Numeric vector. A time series of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior. The code has been optimized to work with integer values. Please re-scale the input values if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000).
#' @param dates A date vector. The number of dates must be equal to the number of "x" values (numeric input vector).
#' @param h Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h = 1 if the vegetation is in the Northern Hemisphere (season starting on January 1st), h = 2 if it is in the Southern Hemisphere (season starting on July 1st).
#' @param xlab Character vector (or expression) giving plot title in x axis label (e.g. xlab = "day of the growing season").
#' @param ylab Character vector (or expression) giving plot title in y axis label (e.g. ylab = "NDVI").
#' @param rge  Numeric vector with the minimum and maximum values of the vegetation index (e.g. NDVI) used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000).
#' @details It is a "heatmap" of the annual phenological variability of the \code{\link{Phen}} output. It calculates and plot a likelihood map of the vegetation-index–time space using a numeric vector of greenness proxies such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). Also a vector with dates for the vegetation index values is required. This function is partially based on the ci2d function on package \href{https://CRAN.R-project.org/package=gplots}{gplots}.
#' @seealso \code{\link{Phen}}, \code{\link{PhenMap}}
#' @examples
#' \dontshow{
#' ## Testing the function with an NDVI time series of a deciduous Nothofagus macrocarpa forest
#' # Load data
#' data("phents")
#' # PhenKplot for the given data
#' PhenKplot(
#'   x = phents$NDVI, dates = phents$dates, h = 2,
#'   xlab = "Day of the growing season",
#'   ylab = "NDVI", rge = c(0, 10000)
#' )
#' }
#' \donttest{
#' library(lubridate)
#' library(terra)
#' ## Testing raster data from Central Chile (NDVI), h=2##
#' # Load data
#' f <- system.file("extdata/MegaDrought_spatRast.rda", package = "npphen")
#' MegaDrought <- readRDS(f)
#' # Dates
#' data("modis_dates")
#'
#' # Generate a Raster time series from a particular pixel 
#' # using a SpatRaster and a date for Central Chile
#' md_pixel <- cellFromXY(MegaDrought, cbind(313395, 6356610))
#' md_pixelts <- as.numeric(MegaDrought[md_pixel])
#' plot(modis_dates, md_pixelts, type = "l")
#'
#' # Variability of the annual phenology for the given pixel
#' PhenKplot(x = md_pixelts, dates = modis_dates,
#' h = 2, xlab = "DGS", ylab = "NDVI", rge = c(0, 10000))
#'
#'
#' ## Testing with the Bdesert_spatRast from 
#' ## the Atacama Desert, Northern Chile (NDVI), h=2 ##
#'
#' # Load data
#' # SparRaster
#' f <- system.file("extdata/Bdesert_spatRast.rda", package = "npphen")
#' Bdesert <- readRDS(f)
#'
#' # Generate a Raster time series from a particular pixel 
#' # using a SpatRaster and a date for Northern Chile
#' bd_pixel <- cellFromXY(Bdesert, cbind(286638, 6852107))
#' bd_pixelts <- as.numeric(Bdesert[bd_pixel])
#' plot(modis_dates, bd_pixelts, type = "l")
#'
#' # Variability of the annual phenology for the given pixel
#' PhenKplot(x = bd_pixelts, dates = modis_dates, 
#' h = 2, xlab = "DGS", ylab = "NDVI", rge = c(0, 10000))
#' }
#' @export

PhenKplot <-
  function(x, dates, h, xlab, ylab, rge) {
    if (length(rge) != 2) {
      stop("rge must be a vector of length 2")
    }
    if (rge[1] > rge[2]) {
      stop("rge vector order must be minimum/maximum")
    }
    if (length(dates) != length(x)) {
      stop("N of dates and files do not match")
    }

    nGS <- 365L

    if (all(is.na(x))) {
      stop("Vector with only NA's. Please check your input data")
    }

    if (all(x < rge[1]) | all(x > rge[2], na.rm = T)) {
      stop("Inconsistency between rge and x. Please check your input data")
    }

    DOY <- lubridate::yday(dates)
    DOY[which(DOY == 366)] <- 365L
    D1 <- cbind(DOY, x)
    if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * nrow(D1))) {
      return(rep(NA, nGS))
    }

    if (h != 1 && h != 2) {
      stop("Invalid h")
    }
    DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
    if (h == 2) {
      for (i in 1:nrow(D1)) {
        D1[, 1] <- DOGS[match(D1[, 1], DOGS[, 1]), 2]
      }
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
    
    first.no.NA.DOY <- min(D1[,1][which(is.na(D1[,2])==FALSE)])
    last.no.NA.DOY <- max(D1[,1][which(is.na(D1[,2])==FALSE)])

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

    h2d <- list()
    h2d$x <- 1L:365L
    h2d$y <- seq(rge[1], rge[2], length.out = 500)
    total_density <- sum(K1Con)
    h2d$density <- K1Con / total_density
    uniqueVals <- rev(unique(sort(h2d$density)))
    cumRFDs <- cumsum(uniqueVals)
    density_indices <- match(h2d$density, uniqueVals)
    h2d$cumDensity <- matrix(cumRFDs[density_indices], nrow = nrow(h2d$density), ncol = ncol(h2d$density))
    
    na.sta <- first.no.NA.DOY-1
    na.end <- last.no.NA.DOY+1
    if(na.sta>=1) {h2d$cumDensity[1:na.sta,] <- NA}
    if(na.end<=365) {h2d$cumDensity[na.end:365,] <- NA}

    image(h2d$x, h2d$y, h2d$cumDensity, xlab = xlab, ylab = ylab, font.lab = 2, breaks = c(0, 0.5, 0.75, 0.9, 0.95), col = grDevices::heat.colors(n = 4, alpha = 0.6))
    contour(h2d$x, h2d$y, h2d$cumDensity, levels = c(0, 0.5, 0.75, 0.9, 0.95), add = T, col = grDevices::grey(0.25), labcex = 1)
    lines(seq(1, 365), MAXY, lwd = 3, col = "dark red")
  }
