# npphen <img src="man/figures/npphen_logo.png" align="right" height=150, width = 150, alt="100" />

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/npphen)](https://CRAN.R-project.org/package=npphen)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/npphen)](https://www.pucv.cl/uuaa/site/edic/base/port/labgrs.html)


<!-- badges: end -->

The R package npphen has been developed primarily to enable land surface phenology reconstruction and anomaly detection by using remote sensing time series. Nevertheless, npphen includes functions to analyze any type of numerical time series. Examples of remote sensing time series that can be tackled using npphen are NDVI (Normalized Difference Vegetation Index) time series from the [GIMMS project](https://glam1.gsfc.nasa.gov/), from the [U.S. MODIS Program](https://modis.gsfc.nasa.gov/data/dataprod/mod13.php) or the E.S.A. [Sentinel 2 Program](https://scihub.copernicus.eu/dhus/#/home). Examples of numerical time series can be the NPP (Net Primary Productivity) from a [FluxNet Tower](https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=9).

## Features

* Long term phenological reconstruction for satellite time series, with no need of gap filling or any other adjustment.
  * Simple phenology
  * Phenology and confidence intervals
* Custom anomaly calculation using the long term phenology as baseline
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

For a brief introduction and examples, check out the [tutorial](https://www.pucv.cl/uuaa/labgrs/proyectos/introduction-to-npphen-in-r). See changelogs on [NEWS.md](https://github.com/labGRS/npphen/blob/master/NEWS.md)
  

## License

The **npphen** package as a whole is licensed under the GPLv3. See the [LICENSE](LICENSE) file for more details.

## References

See the recently published paper for applications and examples use of **npphen**: 

* Chávez, R.O.; Estay, S.A.; Lastra, J.A.; Riquelme, C.G.; Olea, M.; Aguayo, J.; Decuyper, M. **npphen: An R-Package for Detecting and Mapping Extreme Vegetation Anomalies Based on Remotely Sensed Phenological Variability**. Remote Sens. 2023, 15, 73. https://doi.org/10.3390/rs15010073 

* Estay, S.A.; Chávez, R.O.; Rocco, R.; Gutiérrez, A.G. Quantifying Massive Outbreaks of the Defoliator Moth Ormiscodes Amphimone in Deciduous Nothofagus-Dominated Southern Forests Using Remote Sensing Time Series Analysis. J. Appl. Entomol. 2019, 143, 787–796. https://doi.org/10.1111/jen.12643

* Chávez, R.O.; Rocco, R.; Gutiérrez, Á.G.; Dörner, M.; Estay, S.A. A Self-Calibrated Non-Parametric Time Series Analysis Approach for Assessing Insect Defoliation of Broad-Leaved Deciduous Nothofagus Pumilio Forests. Remote Sens. 2019, 11, 204. https://doi.org/10.3390/rs11020204

* Decuyper, M.; Chávez, R.O.; Čufar, K.; Estay, S.A.; Clevers, J.G.P.W.; Prislan, P.; Gričar, J.; Črepinšek, Z.; Merela, M.; de Luis, M.; et al. Spatio-Temporal Assessment of Beech Growth in Relation to Climate Extremes in Slovenia—An Integrated Approach Using Remote Sensing and Tree-Ring Data. Agric. For. Meteorol. 2020, 287, 107925. https://doi.org/10.1016/j.agrformet.2020.107925

* Chávez, R.O.; Moreira-Muñoz, A.; Galleguillos, M.; Olea, M.; Aguayo, J.; Latín, A.; Aguilera-Betti, I.; Muñoz, A.A.; Manríquez, H. GIMMS NDVI Time Series Reveal the Extent, Duration, and Intensity of “Blooming Desert” Events in the Hyper-Arid Atacama Desert, Northern Chile. Int. J. Appl. Earth Obs. Geoinf. 2019, 76, 193–203. https://doi.org/10.1016/j.jag.2018.11.013

* Chávez, R.O.; Castillo-Soto, M.E.; Traipe, K.; Olea, M.; Lastra, J.A.; Quiñones, T. A Probabilistic Multi-Source Remote Sensing Approach to Evaluate Extreme Precursory Drought Conditions of a Wildfire Event in Central Chile. Front. Environ. Sci. 2022, 10, 427. https://doi.org/10.3389/fenvs.2022.865406


