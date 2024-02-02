#' Smooth EMOS (SEMOS)
#'
#' Smooth EMOS (SEMOS) as described by Jobst, Möller and Groß (2024).
#'
#' @param train data frame containing the training data.
#' @param test data frame containing the testing data.
#' @param doy_col column of the variable day of the year.
#' @param obs_col column of the observation variable.
#' @param mean_col column of the ensemble mean forecast variable.
#' @param sd_col column of the ensemble standard deviation variable.
#' @param ... unused
#'
#' @return A data frame containing the distribution (location and scale) parameters.
#'
#' @references
#'
#' Jobst, D., Möller, A., and Groß, J. (2024). Time Series based Ensemble Model Output Statistics for Temperature Forecasts Postprocessing. https://doi.org/10.48550/arXiv.2402.00555.
#'
#' @examples
#' # load data for station Hannover
#' data(station)
#'
#' # select data for lead time 24 hours
#' data <- station[station$lt == 24, ]
#'
#' # split data in training and test data
#' train <- data[data$date <= as.Date("2019-12-31"), ]
#' test <- data[data$date > as.Date("2019-12-31"), ]
#'
#' # estimate model
#' fit <- semos(train = train,
#'              test = test,
#'              doy_col = 3,
#'              obs_col = 9,
#'              mean_col = 10,
#'              sd_col = 11)
#'
#' # distribution parameters
#' head(fit)
#'
#' # select data for lead time 48 hours
#' data <- station[station$lt == 48, ]
#'
#' # split data in training and test data
#' train <- data[data$date <= as.Date("2019-12-31"), ]
#' test <- data[data$date > as.Date("2019-12-31"), ]
#'
#' # estimate model
#' fit <- semos(train = train,
#'              test = test,
#'              doy_col = 3,
#'              obs_col = 9,
#'              mean_col = 10,
#'              sd_col = 11)
#'
#' # distribution parameters
#' head(fit)
#'
#' @importFrom stats ar dnorm pnorm embed lm na.omit optim
#' @importFrom imputeTS na_ma
#' @export
semos <- function(train,
                  test,
                  doy_col,
                  obs_col,
                  mean_col,
                  sd_col,
                  ...) {

  # get necessary training data
  n <- nrow(train)
  doy <- train[, doy_col]
  obs <- train[, obs_col]
  # impute missing values due to leadtime to be able to calculate approximate residuals
  obs <- na_ma(obs)
  m <- train[, mean_col]
  s <- train[, sd_col]

  if(! all((1:365 %in% unique(doy))) ) {
    stop("Each day of the year needs to be at least once in doy!")
  }

  sin1 <- sin(2*pi*doy/365.25)
  sin2 <- sin(4*pi*doy/365.25)
  cos1 <- cos(2*pi*doy/365.25)
  cos2 <- cos(4*pi*doy/365.25)
  doys <- cbind(1, sin1, sin2, cos1, cos2)
  # data frame for seasonal linear model
  df <- data.frame(obs, m, sin1, sin2, cos1, cos2)

  # fit seasonal linear model
  slm <- lm(obs ~ sin1 + sin2 + cos1 + cos2 + (1 + sin1 + sin2 + cos1 + cos2)*m, data = df)

  optim_fun <- function(pars, obs, m, s, doys) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)

      z <- (obs - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
  res <- optim(par = as.numeric(c(slm$coefficients, rep(0, 5), 1, rep(0, 4))),
               fn = optim_fun,
               obs = obs,
               m = m,
               s = s,
               doys = doys,
               method = "BFGS")

    ##############################
    # START PREDICTION
    ##############################

    # get parameter estimates
    pars <- res$par

    # get necessary test data
    doy <- test[, doy_col]
    obs <- test[, obs_col]
    m <- test[, mean_col]
    s <- test[, sd_col]

    sin1 <- sin(2*pi*doy/365.25)
    sin2 <- sin(4*pi*doy/365.25)
    cos1 <- cos(2*pi*doy/365.25)
    cos2 <- cos(4*pi*doy/365.25)
    doys <- cbind(1, sin1, sin2, cos1, cos2)

    # predict values
    MU <- as.vector(doys %*% pars[1:5] + doys %*% pars[6:10]*m)
    SIGMA <- as.vector(exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s))


  # create output
  output <- data.frame(location = MU, scale = SIGMA)

  # output
  return(output)

}
