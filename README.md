
# tsEMOS: Time Series based Ensemble Model Output Statistics

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/tsEMOS)](https://CRAN.R-project.org/package=tsEMOS)
[![R-CMD-check](https://github.com/jobstdavid/tsEMOS/workflows/R-CMD-check/badge.svg)](https://github.com/jobstdavid/tsEMOS/actions)
[![version](https://img.shields.io/badge/version-0.2.0-green.svg?style=flat)](https://github.com/jobstdavid/tsEMOS)
<!-- badges: end -->

An R package for time series based extensions of Ensemble Model Output
Statistics (EMOS) as described in the references.

It depends on the R-packages:

- [imputeTS](https://cran.r-project.org/web/packages/imputeTS/index.html):
  Time series missing value imputation.
- [rugarch](https://cran.r-project.org/web/packages/rugarch/index.html):
  Univariate GARCH models.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("jobstdavid/tsEMOS")
```

## Package overview

Below is an overview of all functions contained in the R-package for
model estimation and prediction:

- `semos`: smooth EMOS (SEMOS).
- `dar_semos`: deseasonalized autoregressive smooth EMOS (DAR-SEMOS).
- `dargarchmult_semos`: multiplicative deseasonalized autoregressive
  smooth EMOS with generalized autoregressive conditional
  heteroscedasticity (DAR-GARCH-SEMOS ($\cdot$)).
- `dargarchadd_semos`: additive deseasonalized autoregressive smooth
  EMOS with generalized autoregressive conditional heteroscedasticity
  (DAR-GARCH-SEMOS (+)).
- `sar_semos`: standardized autoregressive smooth EMOS (SAR-SEMOS).

## Example

### Load R-package and data

``` r
# load package
library(tsEMOS)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo

# load data for station Hannover
data(station)

# select data for lead time 24 hours
data <- station[station$lt == 24, ]

# split data in training and test data
train <- data[data$date <= as.Date("2019-12-31"), ]
test <- data[data$date > as.Date("2019-12-31"), ]
```

### SEMOS

``` r
fit <- semos(train = train,
             test = test,
             doy_col = 3,
             obs_col = 9,
             mean_col = 10,
             sd_col = 11,
             n_ahead = 0)
```

### DAR-SEMOS

``` r
fit <- dar_semos(train = train,
                 test = test,
                 doy_col = 3,
                 obs_col = 9,
                 mean_col = 10,
                 sd_col = 11,
                 n_ahead = 0)
```

### DAR-GARCH-SEMOS ($\cdot$)

``` r
fit <- dargarchmult_semos(train = train,
                          test = test,
                          doy_col = 3,
                          obs_col = 9,
                          mean_col = 10,
                          sd_col = 11,
                          n_ahead = 0)
```

### DAR-GARCH-SEMOS (+)

``` r
fit <- dargarchadd_semos(train = train,
                         test = test,
                         doy_col = 3,
                         obs_col = 9,
                         mean_col = 10,
                         sd_col = 11,
                         n_ahead = 0)
```

### SAR-SEMOS

``` r
fit <- sar_semos(train = train,
                 test = test,
                 doy_col = 3,
                 obs_col = 9,
                 mean_col = 10,
                 sd_col = 11,
                 n_ahead = 0)
```

## Contact

Feel free to contact <jobstd@uni-hildesheim.de> if you have any
questions or suggestions.

## References

Jobst, D., Möller, A., and Groß, J. (2024). Time Series based Ensemble
Model Output Statistics for Temperature Forecasts Postprocessing.
<https://doi.org/10.48550/arXiv.2402.00555>.
