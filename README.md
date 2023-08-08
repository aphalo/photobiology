
# photobiology <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version-last-release/photobiology)](https://cran.r-project.org/package=photobiology)
[![cran
checks](https://badges.cranchecks.info/worst/photobiology.svg)](https://cran.r-project.org/web/checks/check_results_photobiology.html)
[![R build
status](https://github.com/aphalo/photobiology/workflows/R-CMD-check/badge.svg)](https://github.com/aphalo/photobiology/actions)
<!-- badges: end -->

Package ‘**photobiology**’ defines a system of classes for storing
spectral data and accompanying metadata. For each of these classes
specialised summary methods, maths operators and functions are provided.
In addition, classes for storing collections of objects of the classes
for individual spectra are defined as well as ‘apply’ functions.
Extraction and replacement operators are implemented.

Functions for calculation of the position of the sun, times of sunrise
and sunset, day length and night length are also implemented.

The package supports storage and manipulation of data for radiation
quantities and for optical properties of objects and solutes, as well as
response and action spectra for photobiological, photochemical and
photoelectrical responses.

This package is the core of a suite of R packages for photobiological
calculations described at the
[r4photobiology](https://www.r4photobiology.info) web site and in the
vignette [The R for Photobiology
Suite](https://docs.r4photobiology.info/photobiology/articles/userguide-0-r4p-introduction.html).

## Examples

The first example shows you how to estimate solar irradiance in W/m2
under a filter. We use a measured solar spectrum and a measured filter
transmission spectrum.

``` r
library(photobiology)
e_irrad(sun.spct * yellow_gel.spct)
#>  E_Total 
#> 146.5053 
#> attr(,"time.unit")
#> [1] "second"
#> attr(,"radiation.unit")
#> [1] "total energy irradiance"
```

The second example shows some simple astronomical calculations for the
sun.

``` r
geocode <- data.frame(lon = 0, lat = 55)
date <- lubridate::now(tzone = "UTC")
sunrise_time(date, tz = "UTC", geocode = geocode)
#> [1] "2023-08-08 04:21:26 UTC"
day_length(date, tz = "UTC", geocode = geocode)
#> [1] 15.47508
```

## Installation

Installation of the most recent stable version from CRAN:

``` r
install.packages("photobiology")
```

Once package ‘photobiology’ is installed, installation of the remaining
or missing packages in the suite from CRAN:

``` r
intalled_pkgs <- installed.packages()[ , 1]
missing_pkgs <- setdiff(photobiology::r4p_pkgs, intalled_pkgs)
if (length(missing_pkgs) > 0) {
 install.packages(missing_pkgs)
}
```

Installation of the current unstable version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("aphalo/photobiology")
```

## Documentation

HTML documentation is available at
(<https://docs.r4photobiology.info/photobiology/>), including three
*User Guides*.

News on updates to the different packages of the ‘r4photobiology’ suite
are regularly posted at (<https://www.r4photobiology.info/>).

Two articles introduce the basic ideas behind the design of the suite
and describe its use: Aphalo P. J. (2015)
(<https://doi.org/10.19232/uv4pb.2015.1.14>) and Aphalo P. J. (2016)
(<https://doi.org/10.19232/uv4pb.2016.1.15>).

A book is under preparation, and the draft is currently available at
(<https://leanpub.com/r4photobiology/>).

A handbook written before the suite was developed contains useful
information on the quantification and manipulation of ultraviolet and
visible radiation: Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A.
R., Robson, T. M., & Rosenqvist, E. (Eds.) (2012) Beyond the Visible: A
handbook of best practice in plant UV photobiology (1st ed., p. xxx +
174). Helsinki: University of Helsinki, Department of Biosciences,
Division of Plant Biology. ISBN 978-952-10-8363-1 (PDF),
978-952-10-8362-4 (paperback). PDF file available from
(<https://hdl.handle.net/10138/37558>).

## Contributing

Pull requests, bug reports, and feature requests are welcome at
(<https://github.com/aphalo/photobiology>).

## Citation

If you use this package to produce scientific or commercial
publications, please cite according to:

``` r
citation("photobiology")
#> To cite package ‘photobiology’ in publications use:
#> 
#>   Aphalo, Pedro J. (2015) The r4photobiology suite. UV4Plants Bulletin,
#>   2015:1, 21-29. DOI:10.19232/uv4pb.2015.1.14
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     author = {Pedro J. Aphalo},
#>     title = {The r4photobiology suite},
#>     journal = {UV4Plants Bulletin},
#>     volume = {2015},
#>     number = {1},
#>     pages = {21-29},
#>     year = {2015},
#>     doi = {10.19232/uv4pb.2015.1.14},
#>   }
```

## License

© 2012-2023 Pedro J. Aphalo (<pedro.aphalo@helsinki.fi>). Released under
the GPL, version 2 or greater. This software carries no warranty of any
kind.
