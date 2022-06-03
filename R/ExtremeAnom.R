#' @title ExtremeAnom
#' @description Based on the annual reference frequency distribution of a vegetation index time series (e.g. a numeric vector of NDVI), calculates anomalies and how extreme these anomalies are (rfd position ranging from 0 to 100)
#' @encoding UTF-8
#' @param x Numeric vector. A time series of a vegetation index (e.g. LAI, NDVI, EVI) or any other variable with seasonal behavior. The code has been optimized to work with integer values. Please re-scale the input values if necessary (e.g. NDVI ranging from 0.0000 to 1.0000, multiply by 10,000
#' @param dates A date vector. The number of dates must be equal to the number of "x" values (numeric input vector).
#' @param h Numeric. Indicates the geographic hemisphere. This argument defines the starting date of the growing season. h=1 for the Northern Hemisphere (season starting on January 1st), h=2 for the Southern Hemisphere (season starting on July 1st).
#' @param refp Numeric vector with the correlative number of dates to be used as reference period. For example, refp = c(1:388) for MODIS Vegetation Index 16-days composites MOD13Q1 (18/02/2000 – 18/12/2016)
#' @param anop Numeric vector with the correlative number of dates for the period in which the anomalies and rfd position (how extreme the anomalies are) will be calculated. For example anop = c(389:411) for the complete 2017 of MODIS Vegetation Index 16-days composites MOD13Q1 (01/01/2017 – 19/12/2017). anop y refp can be overlapped.
#' @param rge Numeric vector with two values setting the minimum and maximum values of the response variable (e.g. NDVI) used in the analysis. We suggest the use of theoretically based limits. For example in the case of MODIS NDVI or EVI, it ranges from 0 to 10,000, so rge = c(0,10000)
#' @param output Character string. Defines the output values. 'both' (default) returns both anomalies and rfd position together as a single numeric vector, 'anomalies' returns only anomalies, 'rfd' returns only rfd values (how extreme the anomalies are) and 'clean' returns only extreme anomalies, i.e. anomalies at which a given rfd is overpass (e.g. 0.90). This critical threshold is set by the users using the rfd argument.
#' @param rfd Numeric. This argument only applies when the argument output='clean'. It defines the percentile (from 0 to 0.99) of the reference frequency distribution, for which anomalies are not flagged as extreme anomalies. For example, if 'rfd = 0.90' only anomalies falling outside the '0.90 rfd' (default value) will be flagged as extreme anomalies while the rest will be neglected (NA values). Please notice that 'rfd = 0.90' implies that the 5\% of the most extreme positive and 5\% of the most extreme negative anomalies will be considered.

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
#' anom_rfd <- ExtremeAnom(x = phents$NDVI,dates = all.dates,h = 2,refp = c(1:423),
#' anop = c(1:929),rge = c(0,10000),output = 'both',rfd = 0.90)
#' 
#' selection <- which(anom_rfd[930:1858] > 90)
#' 
#' #basic plot
#' barplot(names=format.Date(all.dates[selection],format='%Y-%m'),anom_rfd[selection],col= ifelse(anom_rfd[selection] < 0,"red","blue"),main ='Anomalies whit rfd > 95%')
#' abline(h=0)
#' 
#' }
#'
#' \donttest{
#' library(lubridate)
#'
#' ## Testing with the MegaDrought_stack from Central Chile (NDVI), 
#' #showing extreme negative anomalies (browning)##
#'
#' # Load data
#' #RasterStack
#' data("MegaDrought_stack")
#' #Dates
#' data("modis_dates")
#'
#' # Extracting a time series from a particular pixel
#' md_pixel <- cellFromXY(MegaDrought_stack,c(313395,6356610))
#' md_pixelts <- as.numeric(MegaDrought_stack[md_pixel])
#' plot(modis_dates,md_pixelts, type='l')
#'
#' # Anomaly detection for the given pixel
#' anomRFD_MD <- ExtremeAnom(x = md_pixelts,dates = modis_dates,
#' h = 2,refp = c(1:423), anop = c(884:929),rfd = 0.9,output = 'both',rge = c(0,10000))
#'
#'#Basic plot
#'
#' selection <- which(anomRFD_MD[47:92] > 90)
#' 
#' barplot(names=format.Date(modis_dates[884:929],format='%Y-%m'),
#' anomRFD_MD[1:46],col= ifelse(anomRFD_MD[1:46] < 0,"red","blue"),
#' main ='Anomalies whit rfd > 0.90')
#' abline(h=0)
#'
#' ## Testing with the Bdesert_stack from the Atacama Desert, Northern Chile (NDVI), 
#' #showing extreme positive anomalies (greening)##
#'
#' # Load data
#' #RasterStack
#' data("Bdesert_stack")
#' #Dates
#' data("modis_dates")
#'
#' # Extracting a time series from a particular pixel
#' bd_pixel<-cellFromXY(Bdesert_stack,c(286638,6852107))
#' bd_pixelts<-as.numeric(Bdesert_stack[bd_pixel])
#' plot(modis_dates,bd_pixelts, type = 'l')
#'
#' # Anomaly detection for the given pixel
#' anomRFD_BD <- ExtremeAnom(x = bd_pixelts,dates = modis_dates,
#' h = 2,refp = c(1:423), anop = c(723:768),rfd = 0.9,output = 'both',rge = c(0,10000))
#' 
#'#Basic plot
#'
#' selection <- which(anomRFD_BD[47:92] > 95)
#' 
#' #basic plot
#' barplot(names=format.Date(modis_dates[723:768],format='%Y-%m'),
#' anomRFD_BD[1:46],col= ifelse(anomRFD_BD[1:46] < 0,"red","blue"),
#' main ='Anomalies whit rfd > 0.95')
#' abline(h=0)
#' }
#' @export

ExtremeAnom <- function(x,dates,h,refp,anop,rge, output = 'both',rfd = 0.90){

    if(length(rge)!=2){stop("rge must be a vector of length 2")}
  if(rge[1]>rge[2]){stop("rge vector order must be minimum/maximum")}
  if(length(dates)!=length(x)){stop("N of dates and files do not match")}
  
  output.method <- match(output, c("both", "anomalies", "rfd","clean"))
  
  if (is.na(output.method)){ 
    stop("Invalid output. Must be 'both', 'anomalies', 'rfd' or 'clean'")}
  
  ref.min <- min(refp)
  ref.max <- max(refp)
  ano.min <- min(anop)
  ano.max <- max(anop)
  ano.len <- ano.max-ano.min+1
  len2 <- 2*ano.len
  
  if(ref.min>=ref.max | ano.min>ano.max){stop("for refp or anop, lower value > upper value")}
  
  if (all(is.na(x)) & output == 'both') {
    return(rep(NA,len2))
  }
  if (all(is.na(x)) & output %in% c('clean','anomalies','rfd')) {
    return(rep(NA,ano.len))
  }
  
  DOY <- yday(dates)
  DOY[which(DOY==366)]<-365
  D1<-cbind(DOY[ref.min:ref.max],x[ref.min:ref.max])
  D2<-cbind(DOY[ano.min:ano.max],x[ano.min:ano.max])
  
  if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*nrow(D1)) & output == 'both') {
    return(rep(NA,len2))
  }
  
  if(length(unique(D1[,2]))<10 | (nrow(D1)-sum(is.na(D1)))<(0.1*nrow(D1)) & output %in% c('clean','anomalies','rfd')) {
    return(rep(NA,ano.len))
  }
  
  if (all(is.na(D2[,2])) & output == 'both') {
    return(rep(NA,len2))
  }
  if (all(is.na(D2[,2])) & output %in% c('clean','anomalies','rfd')) {
    return(rep(NA,ano.len))
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
  
  h2d <- list()
  h2d$x <- seq(1,365)
  h2d$y <- seq(rge[1],rge[2],len=500)
  h2d$density <- K1Con/sum(K1Con)
  uniqueVals <- rev(unique(sort(h2d$density)))
  cumRFDs <- cumsum(uniqueVals)
  names(cumRFDs) <- uniqueVals
  h2d$cumDensity <- matrix(nrow = nrow(h2d$density), ncol = ncol(h2d$density))
  h2d$cumDensity[] <- cumRFDs[as.character(h2d$density)]

  if(h==2){
    for(i in 1:nrow(D2)){
      D2[i,1]<-DOGS[which(DOGS[,1]==D2[i,1],arr.ind=TRUE),2]}} 
  
  Anoma <- rep(NA,ano.len)
  for(i in 1:nrow(D2)){
    Anoma[i]<-as.integer(D2[i,2]-MAXY[D2[i,1]])}
  Anoma <- Anoma[1:ano.len]
  names(Anoma) <- paste('anom',dates[ano.min:ano.max],sep = '_')

  rowAnom<-matrix(NA,nrow=nrow(D2),ncol=500)
  for(i in 1:nrow(D2)){
    rowAnom[i,]<-abs(h2d$y-D2[i,2])}
  rowAnom2<-unlist(apply(rowAnom,1,function(x) {if(all(is.na(x))) {NA} else {which.min(x)}}))
  AnomRFD<-rep(NA,nrow(D2))
  for( i in 1:nrow(D2)){
    AnomRFD[i]<-h2d$cumDensity[D2[i,1],rowAnom2[i]]}
  AnomRFDPerc <- round(100*(AnomRFD))
  names(AnomRFDPerc)<- paste('rfd',dates[ano.min:ano.max],sep = '_')

  rfd <- rfd*100
  if(output == 'both'){
  out_data <- c(Anoma,AnomRFDPerc)}
  if(output == 'anomalies'){
  out_data <- Anoma}
  if(output == 'rfd'){
  out_data <- AnomRFDPerc}
  if(output == 'clean'){
    if (rfd == 0 ){ 
      stop("For output = clean, rfd should be a value between 0 - 0.99")}

  p <- which(AnomRFDPerc >= rfd | is.na(AnomRFDPerc))
  aa <- rep(NA, ano.len)
    for (i in p) {
      aa[i] <- Anoma[i]
    }
    names(aa)<- dates[ano.min:ano.max]
    
  out_data <- aa}

  out_data
  
}
