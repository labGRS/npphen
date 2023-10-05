#' @title ExtremeAnoMap
#' @description Based on the annual reference frequency distribution (rfd) of a vegetation index time series (e.g. a SpatRaster of NDVI), calculates anomalies and how extreme these anomalies are (rfd position ranging from 0 to 100).
#' @encoding UTF-8
#' @param s SpatRaster. A time series of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior.The code has been optimized to work with integer values. Please re-scale the input SpatRaster if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000).
#' @param dates A date vector. The number of dates must be equal to the number of layers of the vegetation index SpatRaster.
#' @param h Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h = 1 if the vegetation is in the Northern Hemisphere (season starting on January 1st), h = 2 if it is in the Southern Hemisphere (season starting on July 1st).
#' @param refp Numeric vector with the correlative number of dates to be used as reference period. For example, refp = c(1:393) for MODIS Vegetation Index 16-days composites (18/02/2000 – 06/06/2017)
#' @param anop Numeric vector with the correlative number of dates for the period in which the anomalies will be calculated. For example refp = c(21:43) for the first complete year for MODIS Vegetation Index 16-days composites (01/01/2001 – 19/12/2001). anop y refp can overlap
#' @param rge Numeric vector with the minimum and maximum values of the vegetation index used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000).
#' @param output Character string. Defines the outputs. output = 'both' (default) returns both the vegetation index anomalies and rfd position together as a single numeric vector, output = 'anomalies' returns only anomalies, output = 'rfd' returns only rfd values (how extreme the anomalies are) and output = 'clean' returns only anomalies at which a given rfd is overpass (e.g. 0.90). This critical threshold is set by the users using the rfd argument.
#' @param rfd Numeric. This argument only applies when the argument output = 'clean'. It defines the percentile (from 0 to 0.99) of the reference frequency distribution, for which anomalies are not flagged as extreme anomalies. For example, if 'rfd = 0.90' only anomalies falling outside the '0.90 rfd' (default value) will be flagged as extreme anomalies while the rest will be neglected (NA values). Please notice that 'rfd = 0.90' implies that the 5\% of the most extreme positive and 5\% of the most extreme negative anomalies will be considered.
#' @param nCluster	Numeric. Number of CPU cores to be used for computational calculations
#' @param outname	Character vector with the output directory path and filename with extension or only the filename and extension if the working directory was set. For example outname = "output_anom.tif". See \code{\link[terra]{writeRaster}}
#' @param datatype	Character. Output data type. For example datatype = "INT2S" for signed integer values. See \code{\link[terra]{writeRaster}}
#' @details Similar to \code{\link{ExtremeAnom}}, it calculates phenological anomalies but using a SpatRaster instead of a numeric vector of vegetation canopy greenness values (e.g. Leaf Area Index, LAI) or satellite based greenness proxies such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). For this purpose, it divides the time series (SpatRaster) of vegetation greenness into 2: the reference period, from which the annual phenological cycle is calculated (same as \code{\link{PhenMap}} function), and the observation period, for which we want to calculate anomalies with respect to the annual phenological cycle. Negative anomalies correspond to observed values lower than the reference and positive anomalies to values higher than the reference. This anomalies can be filtered by the position of the observation within the historical rfd. Users can, for example, set 'rfd = 0.95' to consider only anomalies that outside the 95\% rfd of historical records.
#' @return SpatRaster
#' @seealso \code{\link{ExtremeAnom}}
#' @examples
#' \dontrun{
#' ## DEPENDING ON HARDWARE, THIS PROCESS CAN BE HIGHLY TIME CONSUMING##
#'
#' ## Testing with the MegaDrought_stack from Central Chile (NDVI),
#' # showing extreme negative anomalies (browning)##
#'
#' # Load data
#'
#' # SpatRaster
#' library(terra)
#' f <- system.file("extdata/MegaDrought_spatRast.rda", package = "npphen")
#' MegaDrought <- readRDS(f)
#' # Dates
#' data("modis_dates")
#'
#' # Define the number of cores to be use. In this example we use 1
#' nc1 <- 1
#' ExtremeAnoMap(
#'   s = MegaDrought, dates = modis_dates, h = 2, refp = c(1:423),
#'   anop = c(884:929), rfd = 0.9, output = "both", nCluster = nc1, outname = "anomRFD_MD.tif",
#'   datatype = "INT2S", rge = c(0, 10000)
#' )
#' # map_an1 <- rast("anomRFD_MD.tif") # run only for load anomaly brick
#' # plot(map_an1)
#'
#' ## Testing with the Bdesert_stack from the Atacama Desert, Northern Chile (NDVI),
#' # showing extreme positive anomalies (greening)##
#'
#' # Load data
#' # SpatRaster
#' f <- system.file("extdata/Bdesert_spatRast.rda", package = "npphen")
#' Bdesert <- readRDS(f)
#' # Dates
#' data("modis_dates")
#'
#' # Define the number of cores to be use. In this example we use 1
#' nc1 <- 1
#'
#' ExtremeAnoMap(
#'   s = Bdesert, dates = modis_dates, h = 2, refp = c(1:423),
#'   anop = c(723:768), rfd = 0.9, output = "both", nCluster = nc1, outname = "anomRFD_BD.tif",
#'   datatype = "INT2S", rge = c(0, 10000)
#' )
#' # map_an1 <- rast("anomRFD_BD.tif") # run only for load anomaly brick
#' # plot(map_an1)
#' }
#' @export

ExtremeAnoMap <-
  function(s, dates, h, refp, anop, rge, output = "both", rfd = 0.9, nCluster, outname, datatype) {
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
      if ((all(x < rge[1], na.rm = T) & output %in% c("clean", "anomalies", "rfd")) |
        (all(x > rge[2], na.rm = T) & output %in% c("clean", "anomalies", "rfd"))) {
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
        if (output %in% c("clean", "anomalies", "rfd")) {
          return(rep(NA, ano.len))
        }
      }

      if (all(is.na(D2[, 2])) & output == "both") {
        return(rep(NA, len2))
      }
      if (all(is.na(D2[, 2])) & output %in% c("clean", "anomalies", "rfd")) {
        return(rep(NA, ano.len))
      }

      if (h != 1 && h != 2) {
        stop("Invalid h")
      }
      DOGS <- cbind(seq(1, 365), c(seq(185, 365), seq(1, 184)))
      if (h == 2) {
        for (i in 1:nrow(D1)) {
          D1[i, 1] <- DOGS[which(DOGS[, 1] == D1[i, 1], arr.ind = TRUE), 2]
        }
      }

      Hmat <- ks::Hpi(na.omit(D1))
      Hmat[1, 2] <- Hmat[2, 1]
      K1 <- ks::kde(na.omit(D1), H = Hmat, xmin = c(1, rge[1]), xmax = c(365, rge[2]), gridsize = c(365, 500))
      K1Con <- K1$estimate
      for (j in 1:365) {
        Kdiv <- sum(K1$estimate[j, ])
        ifelse(Kdiv == 0, K1Con[j, ] <- 0, K1Con[j, ] <- K1$estimate[j, ] / sum(K1$estimate[j, ]))
      }
      
      first.no.NA.DOY <- min(D1[,1][which(is.na(D1[,2])==FALSE)])
      last.no.NA.DOY <- max(D1[,1][which(is.na(D1[,2])==FALSE)])
      
      MAXY <- apply(K1Con, 1, max)
      for (i in 1:365) {
        n.select <- which(K1Con[i, ] == MAXY[i], arr.ind = TRUE)
        if (length(n.select) > 1) {
          n <- n.select[1]
          MAXY[i] <- NA
        }
        if (length(n.select) == 1) {
          n <- n.select
          MAXY[i] <- median(K1$eval.points[[2]][n])
        }
        if (i < first.no.NA.DOY) {MAXY[i] <- NA}
        if (i > last.no.NA.DOY) {MAXY[i] <- NA}
      }

      h2d <- list()
      h2d$x <- seq(1, 365)
      h2d$y <- seq(rge[1], rge[2], len = 500)
      h2d$density <- K1Con / sum(K1Con)
      uniqueVals <- rev(unique(sort(h2d$density)))
      cumRFDs <- cumsum(uniqueVals)
      names(cumRFDs) <- uniqueVals
      h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
      h2d$cumDensity[] <- cumRFDs[as.character(h2d$density)]

      na.sta <- first.no.NA.DOY-1
      na.end <- last.no.NA.DOY+1
      h2d$cumDensity[1:na.sta,] <- NA
      h2d$cumDensity[na.end:365,] <- NA
      
      if (h == 2) {
        for (i in 1:nrow(D2)) {
          D2[i, 1] <- DOGS[which(DOGS[, 1] == D2[i, 1], arr.ind = TRUE), 2]
        }
      }

      Anoma <- rep(NA, ano.len)
      for (i in 1:nrow(D2)) {
        Anoma[i] <- D2[i, 2] - MAXY[D2[i, 1]]
      }
      Anoma <- Anoma[1:ano.len]
      names(Anoma) <- paste("anom", dates[ano.min:ano.max], sep = "_")

      rowAnom <- matrix(NA, nrow = nrow(D2), ncol = 500)
      for (i in 1:nrow(D2)) {
        rowAnom[i, ] <- abs(h2d$y - D2[i, 2])
      }
      rowAnom2 <- unlist(apply(rowAnom, 1, function(x) {
        if (all(is.na(x))) {
          NA
        } else {
          which.min(x)
        }
      }))
      AnomRFD <- rep(NA, nrow(D2))
      for (i in 1:nrow(D2)) {
        AnomRFD[i] <- h2d$cumDensity[D2[i, 1], rowAnom2[i]]
      }
      AnomRFDPerc <- round(100 * (AnomRFD))
      names(AnomRFDPerc) <- paste("rfd", dates[ano.min:ano.max], sep = "_")

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
    app(s, fun = ff, filename = outname, cores = nCluster, overwrite = T, wopt = list(datatype = datatype))
  }
