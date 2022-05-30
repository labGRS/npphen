# npphen <img src="man/figures/npphen_logo.png" align="right" height=139 alt="" />

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/npphen)](https://CRAN.R-project.org/package=npphen)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/npphen)](https://www.pucv.cl/uuaa/site/edic/base/port/labgrs.html)


<!-- badges: end -->

The R package npphen has been developed primarily to enable land surface phenology reconstruction and anomaly detection by using remote sensing time series. Nevertheless, npphen includes functions to analyze any type of numerical time series. Examples of remote sensing time series that can be tackled using npphen are NDVI (Normalized Difference Vegetation Index) time series from the [GIMMS project](https://glam1.gsfc.nasa.gov/), from the [U.S. MODIS Program](https://modis.gsfc.nasa.gov/data/dataprod/mod13.php) or the E.S.A. [Sentinel 2 Program](https://scihub.copernicus.eu/dhus/#/home). Examples of numerical time series can be the NPP (Net Primary Productivity) from a [FluxNet Tower](https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=9).

## Features

* Long term phenological reconstruction for satellite time series, with no need of gap filling or any other adjustment.
  * Simple phenology
  * Phenology and confidence intervals
* Custom anomly calculation using the long term phenology as baseline
* Raster functions for raster time series
  

## Installation

To install the stable version from CRAN:

```r
install.packages("npphen")
library(npphen)

```

For developers version:


```r
remotes::install_github('labGRS/npphen')
library(npphen)

```

For a brief introduction and examples, check out the [tutorial](https://www.pucv.cl/uuaa/labgrs/proyectos/introduction-to-npphen-in-r).

## License

The **npphen** package as a whole is licensed under the GPLv3. See the [LICENSE](LICENSE) file for more details.
