npphen 2.1.1
===========

* Vectorize portions of the code for speed improve
* Fix minor bugs

npphen 2.0.0
===========

* Remove raster & rgdal dependencies
* Migrate to terra for raster data processing
* All functions now remove data outside the boundaries of the series to avoid extrapolation
* Fix minor bugs

npphen 1.5.2 
===========

* Fix minor bugs

npphen 1.5.1
===========

* Fix NA problems in **ExtremeAnom** and **ExtremeAnoMap**

npphen 1.5.0
===========

* Add two new functions for anomaly detection:
  * **ExtremeAnom**: calculates anomalies and how extreme these anomalies are (rfd position ranging from 0 to 100) for a numeric vector.
  * **ExtremeAnoMap**: calculates anomalies and how extreme these anomalies are (rfd position ranging from 0 to 100) for a raster stack.
* Fix minor bugs on rfd calculations
* update and replace example data for the package

npphen 1.1.1
===========

* Fix warnings and notes from CRAN

npphen 1.1-0
===========
  
* First version of the package on CRAN
