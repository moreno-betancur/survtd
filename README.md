
<!-- README.md is generated from README.Rmd. Please edit that file -->
survtd
======

[![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)

**survtd** is an R package to fit Cox and additive hazards regression models with time-dependent covariates.

**Note:** Please note that the version available on GitHub is the most up-to-date *development* version of the package. A stable version of the package will be available from CRAN once it is released.

Getting Started
---------------

### Prerequisites

The **survtd** package requires the following packages: **lme4**, **survival**, **timereg**, **mice**, **boot**, **ipw**, **stringr**, **MASS**, **Matrix**.

### Installation

The **survtd** package can be installed directly from GitHub using the **devtools** package. To do this you should first check you have devtools installed by executing the following commands from within your R session:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
```

Then execute the following commands to install **survtd**:

``` r
library(devtools)
install_github("moreno-betancur/survtd")
```

### Reference manual

The full reference manual of the package can be dowloaded [here](https://rawgit.com/moreno-betancur/Reference_manuals/master/survtd.pdf).

Example
-------

An example of usage of the package will be added soon.

Bug Reports
-----------

If you find any bugs, please report them via email to [Margarita Moreno-Betancur](mailto:margarita.moreno@mcri.edu.au).

References
----------

Moreno-Betancur M, Carlin JB, Brilleman SL, Tanamas S, Peeters A, Wolfe R. Survival analysis with time-dependent covariates subject to measurement error and missing data: Two-stage joint model using multiple imputation. 2016 (submitted).
