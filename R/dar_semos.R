#' Deseasonalized Autoregressive smooth EMOS (DAR-SEMOS)
#'
#' Deseasonalized autoregressive smooth EMOS (DAR-SEMOS) as described by Jobst, Möller and Groß (2024).
#'
#' @param train data frame containing the training data.
#' @param test data frame containing the testing data.
#' @param doy_col column of the variable day of the year.
#' @param obs_col column of the observation variable.
#' @param mean_col column of the ensemble mean forecast variable.
#' @param sd_col column of the ensemble standard deviation variable.
#' @param n_ahead integer corresponding to the forecast ahead time
#' (0 for ahead times not greater than 24 hours,
#' 1 for ahead times greater than 24 hours and not greater than 48 hours, and so on).
#' @param aic logical; if \code{TRUE} (default), then the Akaike Information Criterion
#' is used to choose the order of the autoregressive model.
#' If \code{FALSE}, the model of order \code{order.max} is fitted.
#' @param order.max maximum order (or order) of model to fit.
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
#' fit <- dar_semos(train = train,
#'                  test = test,
#'                  doy_col = 3,
#'                  obs_col = 9,
#'                  mean_col = 10,
#'                  sd_col = 11,
#'                  n_ahead = 0)
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
#' fit <- dar_semos(train = train,
#'                  test = test,
#'                  doy_col = 3,
#'                  obs_col = 9,
#'                  mean_col = 10,
#'                  sd_col = 11,
#'                  n_ahead = 1)
#'
#' # distribution parameters
#' head(fit)
#'
#' @importFrom stats ar dnorm pnorm embed lm na.omit optim
#' @importFrom imputeTS na_ma
#' @export
dar_semos <- function(train,
                      test,
                      doy_col,
                      obs_col,
                      mean_col,
                      sd_col,
                      n_ahead = 0,
                      aic = TRUE,
                      order.max = NULL,
                      ...) {


  # get necessary training data
  n <- nrow(train)
  doy <- train[, doy_col]
  # replace values which can not be used for n_ahead by NA
  obs <- c(train[1:(n-n_ahead), obs_col], rep(NA, n_ahead))
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

  # in-sample non-standardized residuals
  r <- as.numeric(slm$residuals)

  # fit AR(p)-process on non-standardized residuals
  arm <- ar(r, aic = aic, order.max = order.max)

  # get parameters of AR(p)-process
  a <- as.numeric(arm$ar)
  mu <- as.numeric(arm$x.mean)
  p <- length(a)
  coef_arm <- c(arm$x.mean, arm$ar)

  # prediction function for the AR(p)-process with n_ahead forecast step
  predict_r <- function(r, mu, a, p, n_ahead) {

    # fill matrix with lags
    r_mat <- na.omit(embed(c(rep(NA, p), r), p))
    # adapt matrix row number based on n_ahead
    r_mat <- matrix(r_mat[-c(nrow(r_mat):(nrow(r_mat)-n_ahead)), ], ncol = p)
    # extend r_mat for prediction values base on n_ahead + 1
    r_mat <- cbind(matrix(NA, ncol = n_ahead + 1, nrow = nrow(r_mat)), r_mat)

    # calculate r-values for correction
    for (k in (n_ahead+1):1) {
      r_mat[, k] <- as.vector(mu + (matrix(r_mat[, (k+1):(k+p)], ncol = p) - mu) %*% a)
    }

    return(r_mat[, 1])

  }

  # distinguish between lag p != 0 and p == 0
  # if p == 0 a seasonal EMOS model without an autoregressive correction is estimated
  if (p != 0) {

    # estimates AR-parameters based on a one step ahead prediction for all forecast steps n_ahead
    # i.e. no in-between step prediction of r-values, which is however needed in practice for n_ahead > 0
    # i.e. one assumes that the observed r-values are already the predicted values
    # however for model prediction, in-between step prediction of r-values needs to be carried out for n_ahead > 0

    optim_fun1 <- function(pars, obs, m, s, doys, p) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)
      r <- as.vector((obs-MU))

      # predict r-values for correction
      r_c <- predict_r(r = r, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = 0)

      # adapt length of parameters based on lag p
      OBS <- obs[-c(1:p)]
      MU <- MU[-c(1:p)]
      SIGMA <- SIGMA[-c(1:p)]

      # correct MU
      MU <- r_c+MU
      # calculate CRPS
      z <- (OBS - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = as.numeric(c(slm$coefficients, rep(0, 5), 1, rep(0, 4), coef_arm)),
                 fn = optim_fun1,
                 obs = obs,
                 m = m,
                 s = s,
                 doys = doys,
                 p = p,
                 method = "BFGS")

    ##############################
    # START PREDICTION
    ##############################

    # get parameter estimates
    pars <- res$par

    # get necessary test data
    doy <- c(doy[(n-p+1-n_ahead):n], test[, doy_col])
    obs <- c(obs[(n-p+1-n_ahead):n], test[, obs_col])
    m <- c(m[(n-p+1-n_ahead):n], test[, mean_col])
    s <- c(s[(n-p+1-n_ahead):n], test[, sd_col])

    sin1 <- sin(2*pi*doy/365.25)
    sin2 <- sin(4*pi*doy/365.25)
    cos1 <- cos(2*pi*doy/365.25)
    cos2 <- cos(4*pi*doy/365.25)
    doys <- cbind(1, sin1, sin2, cos1, cos2)

    # predict values
    MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
    SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)
    r <- as.vector((obs-MU))

    # predict r-values for correction
    r_c <- predict_r(r = r, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = n_ahead)

    # predicted/correct parameters
    MU <- r_c + MU[-c(1:(p+n_ahead))]
    SIGMA <- SIGMA[-c(1:(p+n_ahead))]

  } else {

    optim_fun2 <- function(pars, obs, m, s, doys) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)

      z <- (obs - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = as.numeric(c(slm$coefficients, rep(0, 5), 1, rep(0, 4))),
                 fn = optim_fun2,
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

  }

  # create output
  output <- data.frame(location = MU, scale = SIGMA)

  # output
  return(output)

}
