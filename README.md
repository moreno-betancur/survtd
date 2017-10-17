
<!-- README.md is generated from README.Rmd. Please edit that file -->
survtd
======

[![License](License-GPL-image.svg)](http://www.gnu.org/licenses/gpl-3.0.html)

**survtd** is an R package to fit semi-parametric Cox or additive hazards regression models with time-fixed covariates of any type and multiple continuous time-dependent covariates subject to missing data, measurement error or simply observed in discrete time.

All these issues are handled with the two-stage Multiple Imputation for Joint Modeling (MIJM) approach developed by [Moreno-Betancur et al. (2017)](https://academic.oup.com/biostatistics/article-abstract/4461848/Survival-analysis-with-time-dependent-covariates). Briefly, the MIJM approach uses an adapted multiple imputation by chained equations procedure to impute true values of each marker at each event time.

The package also provides easy-to-use implementations of an unadapted version of the approach (unMIJM); a simple two-stage approach (simple2S), based on single imputation of the continuous markers from linear mixed models; and the last observation carried forward (LOCF) approach for time-dependent covariates of any type, which is the traditional method to incoporate these in hazard models. Data simulation functions also included.

**Note:** Please note that the version available on GitHub is the most up-to-date *development* version of the package. A stable version of the package will be available from CRAN once it is released.

Getting Started
---------------

### Prerequisites

The **survtd** package requires the following packages: **lme4**, **survival**, **timereg**, **mice**, **boot**, **ipw**, **stringr**, **MASS**, **Matrix**.

### Installation

The **survtd** package can be installed directly from GitHub using the **devtools** package, which must be first installed from CRAN by the user to proceed. Then it suffices to execute the following commands to install and load **survtd**:

``` r
library(devtools)
install_github("moreno-betancur/survtd")
library(survtd)
```

### Reference manual

The full reference manual of the package can be dowloaded [here](https://rawgit.com/moreno-betancur/Reference_manuals/master/survtd.pdf).

Example
-------

### Generate example dataset

First we use the `simjm` function to generate an example dataset of n=200 individuals, for whom three moderately correlated time-dependent markers were measured with high missingness and low measurement error, each having a strong effect on the hazard of an event following Cox proportional hazards model. These aspects are all controlled through the arguments of the function (see `?simjm`).

``` r
set.seed(33833)

dat <- simjm(n = 200, surv_model = "Cox", marker_model = "RE", MErr = "Low", 
    Miss = "High", effects = "Strong", corr = "Mod")
```

The variables one would observe in a real study are in the first nine columns, so we will discard the rest which relate to further aspects only available due to knowledge of the data generating model in this simulation context. They are there only so that the true parameter estimates can be easily recovered through the `simjm_benchmark` function (of interest to methodologists).

``` r
dat <- dat[, 1:9]
round(head(dat), 2)
#>    ID   tt event Z1    Z2 tj  Yij_1  Yij_2  Yij_3
#> 1   1 2.43     1  1 54.71  0 138.00 136.07 153.32
#> 2   1 2.43     1  1 54.71  1     NA     NA 148.73
#> 3   1 2.43     1  1 54.71  2     NA     NA     NA
#> 11  2 5.55     1  1 45.67  0 130.36 129.00 139.72
#> 12  2 5.55     1  1 45.67  1     NA     NA     NA
#> 13  2 5.55     1  1 45.67  2 146.65     NA     NA
```

The data is in the long format, having the following columns:

-   The subject identifier: `ID`
-   The time to event or censoring: `tt`
-   The event indicator: `event`
-   A binary time-fixed covariate: `Z1`
-   A continuous time-fixed covariate: `Z2`
-   The visit time: `tj`
-   The value of each time-dependent marker for the given individual and visit time: `Yij_1`, `Yij_2` and `Yij_3`

This is the data format required by the `survtd` function.

### Analyse example dataset

We aim to obtain reliable estimates of the effects of markers `Yij_1`, `Yij_2` and `Yij_3` on the time-to-event, these being defined as the (exponential of) regression coefficients in a Cox model, i.e. hazard ratios. The MIJM approach provides such estimates as it accounts for the missing data and measurement error in the markers.

The approach is applied using a single call to the `survtd` function, as follows:

``` r
fitMIJM <- survtd(Surv(time = tt, event = event) ~ Z1 + Z2 + td(Yij_1) + td(Yij_2) + 
    td(Yij_3), data = dat, id = "ID", visit.time = "tj", model = "Cox", method = "MIJM", 
    M = 5, G = 5, time.trend = as.formula("~x+I(x^2)+I(x^3)"))


# Full results in log-hazard scale:
fitMIJM["Results"]
#> $Results
#>              logHR         SE       CIlow      CIupp      p-value
#> Yij_1  0.008370104 0.01069379 -0.01842104 0.03516124 0.4663082371
#> Yij_2 -0.018458456 0.01778742 -0.06574726 0.02883035 0.3518888008
#> Yij_3  0.027373307 0.01799471 -0.02033261 0.07507922 0.1944972995
#> Z1     0.477369426 0.17550178  0.13142547 0.82331338 0.0070669997
#> Z2     0.063101396 0.01567799  0.03114770 0.09505509 0.0003331663

# Estimates and CIs in hazard ratio scale:
res <- round(exp(fitMIJM["Results"][[1]][, c(1, 3, 4)]), 2)
names(res)[1] <- "HR"
res
#>         HR CIlow CIupp
#> Yij_1 1.01  0.98  1.04
#> Yij_2 0.98  0.94  1.03
#> Yij_3 1.03  0.98  1.08
#> Z1    1.61  1.14  2.28
#> Z2    1.07  1.03  1.10
```

The function uses the usual `Surv` formula syntax, except the time-dependent markers need to be specified by wrapping them with `td()`. We could alternatively fit an additive hazards model (`model="Add"`) or use other estimation approaches (`method="unMIJM`, `method="simple2S"`, `method="LOCF`").

Aside from the number of imputations (`M`) and iterations (`G`) to use in the fitting procedure (relevant to MIJM and unMIJM only), the `time.trend` argument allows specifcation of a function of `x` (in this example a cubic polynomial) that determines how time is modeled in the fixed-effects part of the linear mixed models for the time-dependent markers (relevant to methods MIJM, unMIJM and simple2S).

Bug Reports
-----------

If you find any bugs, please report them via email to [Margarita Moreno-Betancur](mailto:margarita.moreno@mcri.edu.au).

References
----------

Moreno-Betancur M, Carlin JB, Brilleman SL, Tanamas S, Peeters A, Wolfe R (2017). [Survival analysis with time-dependent covariates subject to missing data or measurement error: Multiple Imputation for Joint Modeling (MIJM).](https://academic.oup.com/biostatistics/article-abstract/doi/10.1093/biostatistics/kxx046/4461848/Survival-analysis-with-time-dependent-covariates?redirectedFrom=fulltext) *Biostatistics* \[Epub ahead of print 12 Oct 2017\].

Moreno-Betancur M., Chavance M. (2016) [Sensitivity analysis of incomplete longitudinal data departing from the missing at random assumption: Methodology and application in a clinical trial with drop-outs.](http://journals.sagepub.com/doi/abs/10.1177/0962280213490014) *Statistical methods in medical research*, 25 (4), 1471-1489 \[Epub ahead of print, May 22 2013\].
