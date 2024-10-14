#' @title ExtremeAnom
#' @description Based on the annual reference frequency distribution (rfd) of a vegetation index time series (e.g. a numeric vector of NDVI), calculates anomalies and how extreme these anomalies are (rfd position ranging from 0 to 100)
#' @encoding UTF-8
#' @param x Numeric vector. A time series of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior. The code has been optimized to work with integer values. Please re-scale the input values if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000).
#' @param dates A date vector. The number of dates must be equal to the number of "x" values (numeric input vector).
#' @param h Numeric. Indicates the geographic hemisphere. This argument defines the starting date of the growing season. h = 1 for the Northern Hemisphere (season starting on January 1st), h = 2 for the Southern Hemisphere (season starting on July 1st).
#' @param refp Numeric vector with the correlative number of dates to be used as reference period. For example, refp = c(1:388) for MODIS Vegetation Index 16-days composites MOD13Q1 (18/02/2000 – 18/12/2016)
#' @param anop Numeric vector with the correlative number of dates for the period in which the anomalies and rfd position (how extreme the anomalies are) will be calculated. For example anop = c(389:411) for the complete 2017 of MODIS Vegetation Index 16-days composites MOD13Q1 (01/01/2017 – 19/12/2017). anop & refp can be overlapped.
#' @param rge Numeric vector with two values setting the minimum and maximum values of the response variable (e.g. NDVI) used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000)
#' @param output Character string. Defines the output values. 'both' (default) returns both anomalies and rfd position together as a single numeric vector, 'anomalies' returns only anomalies, 'rfd' returns only rfd values (how extreme the anomalies are) and 'clean' returns only extreme anomalies, i.e. anomalies at which a given rfd is overpass (e.g. 0.90). This critical threshold is set by the users using the rfd argument.
#' @param rfd Numeric. This argument only applies when the argument output = 'clean'. It defines the percentile (from 0 to 0.99) of the reference frequency distribution, for which anomalies are not flagged as extreme anomalies. For example, if 'rfd = 0.90' only anomalies falling outside the '0.90 rfd' (default value) will be flagged as extreme anomalies while the rest will be neglected (NA values). Please notice that 'rfd = 0.90' implies that the 5\% of the most extreme positive and 5\% of the most extreme negative anomalies will be considered.

#' @details Calculates anomalies and a probabilistic measure of how extreme the anomalies are using a numeric vector of vegetation canopy greenness, e.g. satellite based Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). For this purpose, it divides the time series (numeric vector) of vegetation greenness into 2: the reference period, from which the annual phenological reference is calculated (same as \code{\link{Phen}} function), and the observation period, for which we want to calculate anomalies. This anomalies can be filtered by the position of the observation within the historical rfd. Users can, for example, set 'rfd = 0.95' to consider only anomalies that outside the 95\% rfd of historical records.
#'
#' @return numerical vector
#' @seealso \code{\link{ExtremeAnoMap}}
#' @examples
#' \dontshow{
#' ## Testing the function with NDVI time series of deciduous N. macrocarpa forest in Central Chile
#' # Load data
#' data("phents")
#'
#' all.dates <- as.Date(phents$dates)
#'
#' # Reference period: 2000 - 2010 (1:423), anomaly detection period: 2010 onwards (424:929)
#' anom_rfd <- ExtremeAnom(
#'   x = phents$NDVI, dates = all.dates, h = 2, refp = c(1:423),
#'   anop = c(1:929), rge = c(0, 10000), output = "both", rfd = 0.90
#' )
#'
#' selection <- which(anom_rfd[930:1858] > 90)
#'
#' # basic plot
#' barplot(names = format.Date(all.dates[selection], format = "%Y-%m"), anom_rfd[selection], col = ifelse(anom_rfd[selection] < 0, "red", "blue"), main = "Anomalies whit rfd > 95%")
#' abline(h = 0)
#' }
#'
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
#' # Generate a Raster time series using a raster stack and a date database from Central Chile
#' # Obtain data from a particular pixel generating a time series
#' md_pixel <- cellFromXY(MegaDrought, cbind(313395, 6356610))
#' md_pixelts <- as.numeric(MegaDrought[md_pixel])
#' plot(modis_dates, md_pixelts, type = "l")
#'
#' # Anomaly detection for the given pixel
#' anomRFD_MD <- ExtremeAnom(
#'   x = md_pixelts, dates = modis_dates,
#'   h = 2, refp = c(1:423), anop = c(884:929), rfd = 0.9, output = "both", rge = c(0, 10000)
#' )
#'
#' # Basic plot
#'
#' selection <- which(anomRFD_MD[47:92] > 90)
#'
#' barplot(
#'   names = format.Date(modis_dates[884:929], format = "%Y-%m"),
#'   anomRFD_MD[1:46], col = ifelse(anomRFD_MD[1:46] < 0, "red", "blue"),
#'   main = "Anomalies whit rfd > 0.90"
#' )
#' abline(h = 0)
#'
#' ## Testing with the Bdesert_stack from the Atacama Desert, Northern Chile (NDVI),
#' # showing extreme positive anomalies (greening)##
#'
#' # Load data
#' # SparRaster
#' f <- system.file("extdata/Bdesert_spatRast.rda", package = "npphen")
#' Bdesert <- readRDS(f)
#'
#' # Generate a Raster time series using a raster stack and a date database from Northern Chile
#' # Obtain data from a particular pixel generating a time series
#' bd_pixel <- cellFromXY(Bdesert, cbind(286638, 6852107))
#' bd_pixelts <- as.numeric(Bdesert[bd_pixel])
#' plot(modis_dates, bd_pixelts, type = "l")
#'
#' # Anomaly detection for the given pixel
#' anomRFD_BD <- ExtremeAnom(
#'   x = bd_pixelts, dates = modis_dates,
#'   h = 2, refp = c(1:423), anop = c(723:768), rfd = 0.9, output = "both", rge = c(0, 10000)
#' )
#'
#' # Basic plot
#'
#' selection <- which(anomRFD_BD[47:92] > 95)
#'
#' # basic plot
#' barplot(
#'   names = format.Date(modis_dates[723:768], format = "%Y-%m"),
#'   anomRFD_BD[1:46], col = ifelse(anomRFD_BD[1:46] < 0, "red", "blue"),
#'   main = "Anomalies whit rfd > 0.95"
#' )
#' abline(h = 0)
#' }
#' @export

ExtremeAnom <- function(x, dates, h, refp, anop, rge, output = "both", rfd = 0.90) {
  if (length(rge) != 2) {
    stop("rge must be a vector of length 2")
  }
  if (rge[1] > rge[2]) {
    stop("rge vector order must be minimum/maximum")
  }
  if (length(dates) != length(x)) {
    stop("N of dates and files do not match")
  }
  if (length(x) < length(refp) | length(x) < length(anop)) {
    stop("Inconsistent anop or refp. Arguments refp and anop can't be grater than length(x)")
  }
  output.method <- match(output, c("both", "anomalies", "rfd", "clean"))

  if (is.na(output.method)) {
    stop("Invalid output. Must be 'both', 'anomalies', 'rfd' or 'clean'")
  }

  ref.min <- min(refp)
  ref.max <- max(refp)
  ano.min <- min(anop)
  ano.max <- max(anop)
  ano.len <- ano.max - ano.min + 1
  len2 <- 2 * ano.len

  if (ref.min >= ref.max | ano.min > ano.max) {
    stop("for refp or anop, lower value > upper value")
  }

  if (all(is.na(x)) & output == "both") {
    return(rep(NA, len2))
  }
  if (all(is.na(x)) & output %in% c("clean", "anomalies", "rfd")) {
    return(rep(NA, ano.len))
  }

  if ((all(x < rge[1], na.rm = T) & output == "both") | (all(x > rge[2], na.rm = T) & output == "both")) {
    return(rep(NA, len2))
  }
  if ((all(x < rge[1], na.rm = T) & output != 'both') | (all(x > rge[2], na.rm = T) & output != 'both')) {
    return(rep(NA, ano.len))
  }

  DOY <- lubridate::yday(dates)
  DOY[which(DOY == 366)] <- 365
  D1 <- cbind(DOY[ref.min:ref.max], x[ref.min:ref.max])
  D2 <- cbind(DOY[ano.min:ano.max], x[ano.min:ano.max])

  if (length(unique(D1[, 2])) < 10 | (nrow(D1) - sum(is.na(D1))) < (0.1 * nrow(D1))) {
    if (output == "both") {
      return(rep(NA, len2))
    }
    if (output != 'both') {
      return(rep(NA, ano.len))
    }
  }

  if (all(is.na(D2[, 2])) & output == "both") {
    return(rep(NA, len2))
  }
  if (all(is.na(D2[, 2])) & output != 'both') {
    return(rep(NA, ano.len))
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
  }, n.select = n.selects, i = 1:365)

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

  if (h == 2) {
    D2[, 1] <- DOGS[match(D2[, 1], DOGS[, 1]), 2]
  }

  Anoma <- D2[, 2] - MAXY[D2[, 1]]
  names(Anoma) <- paste("anom", dates[ano.min:ano.max], sep = "_")

  rowAnom <- abs(outer(D2[, 2], h2d$y, "-"))
  rowAnom2 <- apply(rowAnom_v2, 1, function(x) {
    if (all(is.na(x))) {
      return(NA) 
    } else {
      return(which.min(x)) 
    }
  })
  AnomRFD <- round(h2d$cumDensity[cbind(D2[, 1], rowAnom2)]*100)
  names(AnomRFD) <- paste("rfd", dates[ano.min:ano.max], sep = "_")

  rfd <- rfd * 100
  if (output == "both") {
    out_data <- c(Anoma, AnomRFDPerc)
  }
  if (output == "anomalies") {
    out_data <- Anoma
  }
  if (output == "rfd") {
    out_data <- AnomRFDPerc
  }
  if (output == "clean") {
    if (rfd == 0) {
      stop("For output = clean, rfd should be a value between 0 - 0.99")
    }

    p <- which(AnomRFDPerc >= rfd | is.na(AnomRFDPerc))
    aa <- rep(NA, ano.len)
    for (i in p) {
      aa[i] <- Anoma[i]
    }
    names(aa) <- dates[ano.min:ano.max]

    out_data <- aa
  }

  out_data
}
