#' Standardized Autoregressive smooth EMOS (SAR-SEMOS)
#'
#' Standardized autoregressive smooth EMOS (SAR-SEMOS) as described by Jobst, Möller and Groß (2024).
#'
#' @param train data frame containing the training data.
#' @param test data frame containing the testing data.
#' @param doy_col column of the variable day of the year.
#' @param obs_col column of the observation variable.
#' @param mean_col column of the ensemble mean forecast variable.
#' @param sd_col column of the ensemble standard deviation variable.
#' @param w integer; window size for the calculation of the empirical standard deviation.
#' Default: \code{w = 1}.
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
#' fit <- sar_semos(train = train,
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
#' fit <- sar_semos(train = train,
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
#' @importFrom stats ar dnorm pnorm embed lm na.omit optim coef predict
#' @importFrom imputeTS na_ma
#' @export
sar_semos <- function(train,
                      test,
                      doy_col,
                      obs_col,
                      mean_col,
                      sd_col,
                      w = 1,
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

  # data frame for seasonal linear model for mu
  df_mu <- data.frame(obs, m, s, sin1, sin2, cos1, cos2)

  # fit seasonal linear model to have initial values for mu
  mu_model <- lm(obs ~ sin1 + sin2 + cos1 + cos2 + (1 + sin1 + sin2 + cos1 + cos2)*m, data = df_mu)

  # calculate observed empirical standard deviation of observation for each doy
  sd.r <- function(r) {
    n <- length(r)
    sqrt(1/(n-1)*sum(r^2))
  }
  # create datafame with residuals
  df_r <- data.frame(id = seq_along(doy), doy = doy, r = mu_model$residuals)

  # extend data.frame with empirical standard deviation of observation for each doy by symmetric window size w
  # enables estimation if only one year of data is available
  s_obs <- c()
  for (k in sort(unique(df_r$doy))) {
    days <- (k + -w:w) %% 365
    days[days == 0] <- 365
    s_obs <- c(s_obs, sd.r(df_r[days, "r"]))
  }
  df_s_obs <- data.frame(doy = sort(unique(df_r$doy)), s_obs = s_obs)

  # combine both data.frame's by doy
  df_s_obs <- merge(df_r, df_s_obs, by = "doy")
  # get correct order of s_obs depending on doy
  df_s_obs <- df_s_obs[order(df_s_obs$id), ]

  # dataframe for seasonal linear model for sigma
  df_sigma <- data.frame(df_mu, s_obs = df_s_obs$s_obs)
  # fit seasonal linear model to have initial values for sigma
  sigma_model <- lm(log(s_obs) ~ sin1 + sin2 + cos1 + cos2 + (1 +  sin1 + sin2 + cos1 + cos2)*s, data = df_sigma)
  # predict the observed standard deviation based on the harmonic model and use link function exp for back-transformation
  s_pred <- exp(predict(sigma_model))

  # get coefficients of linear models mu_model and sigma_model
  coef_lm <- as.numeric(c(coef(mu_model), coef(sigma_model)))

  # in-sample standardized z-residuals
  z <- as.numeric(mu_model$residuals)/s_pred

  # fit AR(p)-process on z-residuals
  arm <- ar(z, aic = aic, order.max = order.max)

  # get parameters of AR(p)-process
  a <- as.numeric(arm$ar)
  mu <- as.numeric(arm$x.mean)
  p <- length(a)
  coef_arm <- as.numeric(c(arm$x.mean, arm$ar))

  # prediction function for the AR(p)-process with n_ahead forecast step
  predict_z <- function(z, mu, a, p, n_ahead) {

    # fill matrix with lags
    z_mat <- na.omit(embed(c(rep(NA, p), z), p))
    # adapt matrix row number based on n_ahead
    z_mat <- matrix(z_mat[-c(nrow(z_mat):(nrow(z_mat)-n_ahead)), ], ncol = p)
    # extend z_mat for prediction values base on n_ahead + 1
    z_mat <- cbind(matrix(NA, ncol = n_ahead + 1, nrow = nrow(z_mat)), z_mat)

    # calculate z-values for correction
    for (k in (n_ahead+1):1) {
      z_mat[, k] <- as.vector(mu + (matrix(z_mat[, (k+1):(k+p)], ncol = p) - mu) %*% a)
    }

    return(z_mat[, 1])

  }

  # distinguish between lag p != 0 and p == 0
  # if p == 0 a seasonal EMOS model without an autoregressive correction is estimated
  if (p != 0) {

    # estimates AR-parameters based on a one step ahead prediction for all forecast steps n_ahead
    # i.e. no in-between step prediction of z-values, which is however needed in practice for n_ahead > 0
    # i.e. one assumes that the observed z-values are already the predicted z-values
    # however for model prediction, in-between step prediction of z-values needs to be carried out for n_ahead > 0

    optim_fun1 <- function(pars, obs, m, s, doys, p) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)
      z <- as.vector((obs-MU)/SIGMA)

      # predict z-values for correction
      z_c <- predict_z(z = z, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = 0)

      # adapt length of parameters based on lag p
      OBS <- obs[-c(1:p)]
      MU <- MU[-c(1:p)]
      SIGMA <- SIGMA[-c(1:p)]

      # correct MU
      MU <- z_c*SIGMA + MU
      # calculate CRPS
      z <- (OBS - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = c(coef_lm, coef_arm),
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
    z <- as.vector((obs-MU)/SIGMA)

    # predict z-values for correction
    z_c <- predict_z(z = z, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = n_ahead)

    # predicted/correct parameters
    MU <- z_c*SIGMA[-c(1:(p+n_ahead))] + MU[-c(1:(p+n_ahead))]
    SIGMA <- SIGMA[-c(1:(p+n_ahead))]

  } else {

    optim_fun2 <- function(pars, obs, m, s, doys) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)

      z <- (obs - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = coef_lm,
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
