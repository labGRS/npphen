#' @title Phen
#' @description Estimates the annual phenological cycle from a time series of vegetation greenness.
#' @encoding UTF-8
#' @param x Numeric vector. A time series of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior. The code has been optimized to work with integer values. Please re-scale the input values if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000
#' @param dates A date vector. The number of dates must be equal to the number of "x" values (numeric input vector).
#' @param h Numeric. Indicates the geographic hemisphere to define the starting date of the growing season. h=1 if the vegetation is in the Northern Hemisphere (season starting at January 1st), h=2 if it is in the Southern Hemisphere (season starting at July 1st)
#' @param frequency Character string. Defines the number of samples for the output phenology and must be one of the this: 'daily' giving output vector of length 365, '8-days' giving output vector of length 46 (i.e MOD13Q1 and MYD13Q1), 'monthly' giving output vector of length 12,'bi-weekly' giving output vector of length 24 (i.e. GIMMS) or '16-days' (default) giving output vector of length 23 (i.e MOD13Q1 or MYD13Q1).
#' @param rge Numeric vector with two values setting the minimum and maximum values of the response variable (e.g. NDVI) used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000)
#' @details Derives the annual phenological cycle for a standard growing season using a numeric vector of vegetation canopy greenness values (e.g. Leaf Area Index, LAI) or satellite based greenness proxies such as the Normalized Difference Vegetation Index (NDVI) or Enhanced Vegetation Index (EVI). A vector with dates for the greenness values is also required.
#' @return A numeric vector, where each value represents the expected greeness at that date
#' @seealso \code{\link{PhenMap}},\code{\link{PhenKplot}}
#' @examples
#' \dontshow{
#' ## Testing function with time series of Nothofagus macrocarpa (NDVI)
#' # Load data
#' data("phents")
#' # Phenology for the given data
#' Phen(x=phents$NDVI,dates=phents$dates,h=2,frequency = '16-days',rge=c(0,10000))
#' }
#' \donttest{
#' library(lubridate)
#'
#' ## Testing raster data from Central Chile (NDVI), h=2##
#'
#' # Load data
#' #RasterStack
#' data("MegaDrought_stack")
#' #Dates
#' data("modis_dates")
#'
#' #Creates Raster time series using a raster stack and a date database from Central Chile
#' # Obtain data from a particular pixel generating a time series
#' md_pixel <- cellFromXY(MegaDrought_stack,c(313395,6356610))
#' md_pixelts <- as.numeric(MegaDrought_stack[md_pixel])
#' plot(modis_dates,md_pixelts, type='l')
#'
#' # Phenology for the given pixel
#' Phen(x=md_pixelts,dates=modis_dates,h=2,frequency='16-days',rge=c(0,10000))
#'
#' ## Testing with the Bdesert_stack from the Atacama Desert, Northern Chile (NDVI), h=2 ##
#'
#' # Load data
#' #RasterStack
#' data("Bdesert_stack")
#' #Dates
#' data("modis_dates")
#' 
#' #Creates Raster time series using a raster stack and a date database from Northern Chile
#' # Obtain data from a particular pixel generating a time series
#' bd_pixel<-cellFromXY(Bdesert_stack,c(286638,6852107))
#' bd_pixelts<-as.numeric(Bdesert_stack[bd_pixel])
#' plot(modis_dates,bd_pixelts, type = 'l')
#'
#' # Phenology for the given pixel
#' Phen(x=bd_pixelts,dates=modis_dates,h=2,frequency='16-days',rge=c(0,10000))
#' }
#' @import raster
#' @import ks
#' @import grDevices
#' @import methods
#' @import rgdal
#' @import snow
#' @importFrom lubridate yday
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom graphics abline axis
#' @export

Phen <-
  function(x,dates,h,frequency = '16-days',rge) {

    if(length(rge)!=2){stop("rge must be a vector of length 2")}
    if(rge[1]>rge[2]){stop("rge vector order must be minimum/maximum")}
    if(length(dates)!=length(x)){stop("N of dates and files do not match")}
    
    freq.method <- match(frequency, c("daily", "8-days", "16-days", 
                               "bi-weekly", "monthly"))
    if (is.na(freq.method)){ 
      stop("Invalid frequency. Must be one of: daily, 8-days, 16-days, bi-weekly or monthly")}
    if(frequency == 'daily'){
      nGS <- 365
    }
    if(frequency == '8-days'){
      nGS <- 46
    }
    if(frequency == '16-days'){
      nGS <- 23
    }
    if(frequency == 'monthly'){
      nGS <- 12
    }
    if(frequency == 'bi-monthly'){
      nGS <- 24
    }
    if (all(is.na(x))) {
      return(rep(NA,nGS))
    }
    
    DOY <- yday(dates)
    DOY[which(DOY==366)]<-365
    D1<-cbind(DOY,x)
    if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*nrow(D1))) {
      return(rep(NA,nGS))
    }

    if(h!=1 && h!=2){stop("Invalid h")}
    DOGS<-cbind(seq(1,365),c(seq(185,365),seq(1,184)))
    if(h==2){
      for(i in 1:nrow(D1)){
        D1[i,1]<-DOGS[which(DOGS[,1]==D1[i,1],arr.ind=TRUE),2]}}
    
    Hmat<-Hpi(na.omit(D1))
    Hmat[1,2]<-Hmat[2,1]
    K1<-kde(na.omit(D1),H=Hmat,xmin=c(1,rge[1]),xmax=c(365,rge[2]),gridsize=c(365,500))
    K1Con<-K1$estimate
    
    for(j in 1:365){ 
      Kdiv<-sum(K1$estimate[j,])
      ifelse(Kdiv==0,K1Con[j,]<-0,K1Con[j,]<-K1$estimate[j,]/sum(K1$estimate[j,]))}
    
    MAXY<-apply(K1Con,1,max)
    for(i in 1:365){ 
      n.select <- which(K1Con[i,]==MAXY[i],arr.ind=TRUE)
      if(length(n.select) > 1){
        n <- n.select[1]
        MAXY[i]<-NA}
      if(length(n.select) == 1){
        n <- n.select
        MAXY[i]<-median(K1$eval.points[[2]][n])
      }
    }

    if(frequency == 'daily'){
      select_DGS <- seq(1,365,1)
      Ref <- MAXY
    }
    if(frequency == '8-days'){
      select_DGS <- seq(1,365,8)
      Ref <- MAXY[select_DGS]
    }
    if(frequency == '16-days'){
      select_DGS <- seq(1,365,16)
      Ref <- MAXY[select_DGS]
    }
    if(frequency == 'monthly'){
      select_DGS <- c(15,46,74,105,135,166,196,227,258,288,319,349)
      Ref <- MAXY[select_DGS]
    }
    if(frequency == 'bi-weekly'){
      select_DGS <- c(1,15,32,46,60,74,91,105,121,135,152,166,182,196,213,227,244,258,274,288,305,319,335,349)
      Ref <- MAXY[select_DGS]
    }
    if(h==1){id.label <- 'DOY'}
    if(h==2){id.label <- 'DGS'}
    plot(select_DGS,Ref,xlab=id.label,ylab='VI',font.lab=2,type='l') 
    axis(1, at = seq(0,365,50))
    
    names(Ref) <- select_DGS 
    Ref
  }

