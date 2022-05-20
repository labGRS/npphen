#' Mega drought in Central Chile
#'
#' MODIS NDVI RasterStack of Central Chile combining MOD13Q1 & MYD13Q1 between 2000.02.18 and 2021.06.26, 
#' cleaned using the detailed QA band. Megadrought affected Central Chile between 2010 and 2015 and again in 2019-2022
#'
#' @format RasterStack EPSG:32719
#' @source \url{https://lpdaac.usgs.gov/}
#' 
"MegaDrought_stack"

#' Blooming Desert Northern Chile
#'
#' MODIS NDVI RasterStack of the Atacama Desert, Northern Chile, combining MOD13Q1 & MYD13Q1 between 2000.02.18 and 2021.06.26,
#' cleaned using the detailed QA band. Blooming deserts are rare greening events to be seen as positive NDVI anomalies in Atacama
#'
#' @format RasterStack EPSG:32719
#' @source \url{https://lpdaac.usgs.gov/}
#' 
"Bdesert_stack"

#' A vector with dates for both RasterStack datasets: Mega drought and Blooming Desert, with 929 dates
#'
#' @format Vector of dates
#' @source \url{https://lpdaac.usgs.gov/}
#' 
"modis_dates"

#' A dataframe with 8-days MODIS NDVI of a pixel with deciduous Nothofagus macrocarpa forest, Central Chile,
#' using MOD13Q1/MYD13Q1 combined 
#'
#'@format A data frame with 929 rows and 2 variables:
  #' \describe{
  #'   \item{dates}{dates for each record}
  #'   \item{NDVI}{NDVI value for a pixel of deciduous N. macrocarpa forest}
  #' }
#' @source \url{https://lpdaac.usgs.gov/}
#' 
"phents"