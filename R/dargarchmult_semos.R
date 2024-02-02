#' Multiplicative Deseasonalized Autoregressive smooth EMOS with Generalized Autoregressive Conditional Heteroscedasticity (DAR-GARCH-SEMOS (*))
#'
#' Multiplicative Deseasonalized Autoregressive smooth EMOS with Generalized Autoregressive
#' Conditional Heteroscedasticity (DAR-GARCH-SEMOS (*)) as described by Jobst, Möller and Groß (2024).
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
#' fit <- dargarchmult_semos(train = train,
#'                           test = test,
#'                           doy_col = 3,
#'                           obs_col = 9,
#'                           mean_col = 10,
#'                           sd_col = 11,
#'                           n_ahead = 0)
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
#' fit <- dargarchmult_semos(train = train,
#'                           test = test,
#'                           doy_col = 3,
#'                           obs_col = 9,
#'                           mean_col = 10,
#'                           sd_col = 11,
#'                           n_ahead = 1)
#'
#' # distribution parameters
#' head(fit)
#'
#' @importFrom stats ar dnorm pnorm embed lm na.omit optim predict residuals
#' @importFrom utils tail
#' @importFrom imputeTS na_ma
#' @importFrom stats4 coef
#' @importFrom rugarch ugarchspec ugarchfit
#'
#' @export
dargarchmult_semos <- function(train,
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
  # impute missing values to be able to calculate approximate residuals of residuals
  obs <- na_ma(obs)
  m <- train[, mean_col]
  s <- train[, sd_col]

  if(! all((1:365 %in% unique(doy))) ) {
    stop("Each day of the year needs to be at least once in doy_col!")
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
  r <- mu_model$residuals

  # calculate observed empirical standard deviation of observation for each doy
  sd.r <- function(r) {
    n <- length(r)
    sqrt(1/(n-1)*sum(r^2))
  }
  # create data.frame with residuals
  df_r <- data.frame(id = seq_along(doy), doy = doy, r = r)

  # extend data.frame with empirical standard deviation of observation for each doy by symmetric window size
  # enables estimation if only one year of data is available
  s_obs <- c()
  for (k in sort(unique(df_r$doy))) {
    days <- (k + -w:w) %% 365
    days[days == 0] <- 365
    s_obs <- c(s_obs, sd.r(df_r[days, "r"]))
  }
  df_s_obs <- data.frame(doy = sort(unique(df_r$doy)), s_obs = s_obs)

  # combine both dataframe's by doy
  df_s_obs <- merge(df_r, df_s_obs, by = "doy")
  # get correct order of s_obs depending on doy
  df_s_obs <- df_s_obs[order(df_s_obs$id), ]

  # data.frame for seasonal linear model for sigma
  df_sigma <- data.frame(df_mu, s_obs = df_s_obs$s_obs)
  # fit seasonal linear model to have initial values for sigma
  sigma_model <- lm(log(s_obs) ~ sin1 + sin2 + cos1 + cos2 + (1 +  sin1 + sin2 + cos1 + cos2)*s, data = df_sigma)
  # predict the observed standard deviation based on the harmonic model and use link function exp for back-transformation
  s_pred <- exp(predict(sigma_model))

  # fit AR(p)-process on non-standardized residuals
  arm <- ar(r, aic = aic, order.max = order.max)

  # get parameters of AR(p)-process
  a <- as.numeric(arm$ar)
  mu <- as.numeric(arm$x.mean)
  p <- length(a)
  coef_arm <- c(arm$x.mean, arm$ar)

  # calculate residuals of autoregressive process
  r_arm <- residuals(arm)
  # calculate deseasonalized residuals
  z <- as.numeric(na.omit(r_arm/s_pred))

  # define GARCH(1,1) model
  spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = (list(armaOrder = c(0,0), include.mean = FALSE)),
                     distribution.model = "norm")

  # fit GARCH(1,1) model
  fit_garch <- ugarchfit(spec, data = z, solver = "hybrid")

  # receive coefficients of GARCH(1,1) model
  omega0 <- as.numeric(coef(fit_garch)[1])
  omega1 <- as.numeric(coef(fit_garch)[3])
  omega2 <- as.numeric(coef(fit_garch)[2])
  coef_garch <- c(omega0, omega1, omega2)

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
  # if p == 0 a seasonal EMOS model without an AR(p) and GARCH(1,1) model is estimated
  if (p != 0) {

    # estimates AR-parameters based on a one step ahead prediction for all forecast steps n_ahead
    # i.e. no in-between step prediction of r-values, which is however needed in practice for n_ahead > 0
    # i.e. one assumes that the observed r-values are already the predicted values
    # however for model prediction, in-between step prediction of r-values needs to be carried out for n_ahead > 0

    optim_fun1 <- function(pars, obs, m, s, doys, p) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- as.vector(exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s))
      r <- as.vector((obs-MU))

      # predict r-values for correction
      r_c <- predict_r(r = r, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = 0)
      n <- length(r_c)
      z <- (tail(r, n = n) - r_c)/tail(SIGMA, n = n)

      # predict z-value for correction
      sig_c <- sqrt(mean(z^2))
      for (k in 1:(n-1)) {
        sig_c <- c(sig_c, sqrt(pars[22+p]^2 + pars[23+p]^2 * sig_c[k]^2 + pars[24+p]^2 * z[k]^2))
      }

      # adapt length of parameters based on lag p
      OBS <- obs[-c(1:p)]
      MU <- MU[-c(1:p)]
      SIGMA <- SIGMA[-c(1:p)]

      # correct MU and SIGMA
      MU <- r_c + MU
      SIGMA <- sig_c * SIGMA
      # calculate CRPS
      z <- (OBS - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = as.numeric(c(mu_model$coefficients, sigma_model$coefficients, coef_arm, coef_garch)),
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

    # get parameters
    pars <- res$par

    # get initial values for sig_c and z from the last iteration in the above optimization
    MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
    SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)
    r <- as.vector((obs-MU))

    r_c <- predict_r(r = r, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = 0)
    n_r_c <- length(r_c)
    z <- (tail(r, n = n_r_c) - r_c)/tail(SIGMA, n = n_r_c)

    sig_c <- sqrt(mean(z^2))
    # predict sig_c
    for (k in 1:(n_r_c-1)) {
      sig_c <- c(sig_c, sqrt(pars[22+p]^2 + pars[23+p]^2 * sig_c[k]^2 + pars[24+p]^2 * z[k]^2))
    }

    # initial sig_c and z
    sig_c <- tail(sig_c, n = 1)
    z <- tail(z, n = 1)

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


    MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
    SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)
    r <- as.vector((obs-MU))
    # predict r-values for correction
    r_c <- predict_r(r = r, mu = pars[21], a = pars[22:(21+p)], p = p, n_ahead = n_ahead)
    n_r_c <- length(r_c)
    z <- c(z, # initial z from above
           (tail(r, n = n_r_c) - r_c)/tail(SIGMA, n = n_r_c))

    # predict sig_c
    for (k in 1:n_r_c) {
      sig_c <- c(sig_c, sqrt(pars[22+p]^2 + pars[23+p]^2 * sig_c[k]^2 + pars[24+p]^2 * z[k]^2))
    }
    # remove initial sig_c
    sig_c <- sig_c[-1]

    # predicted/correct parameters
    MU <- r_c + MU[-c(1:(p+n_ahead))]
    SIGMA <- sig_c * SIGMA[-c(1:(p+n_ahead))]

  } else {

    optim_fun2 <- function(pars, obs, m, s, doys) {

      MU <- doys %*% pars[1:5] + doys %*% pars[6:10]*m
      SIGMA <- exp(doys %*% pars[11:15] + doys %*% pars[16:20]*s)

      z <- (obs - MU)/SIGMA
      crps <- SIGMA * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
      sum(crps)

    }
    res <- optim(par = as.numeric(c(mu_model$coefficients, rep(0, 5), 1, rep(0, 4))),
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

  output <- data.frame(location = MU, scale = SIGMA)

  # output
  return(output)

}
