library("ggfortify")
library("moments")
library(ggplot2)
library(readxl)
library(zoo)
library(xts)
library(gridExtra)
library(rugarch)
library(quantmod)
library(sandwich)
library(lmtest)
library(stringr)
library(readr)
library(dplyr)
library(lubridate)
library(TTR)
library(MSGARCH)
library(tidyr)

set.seed(42) # for filtered hist. sim.

# ------------------------------
# Returns for Model fitting
# ------------------------------

# Load processed data
btc_full <- read.csv("data/btc_full.csv")
btc_full <- btc_full %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))

eth_full <- read.csv("data/eth_full.csv")
eth_full <- eth_full %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))


# Define split date
split_date <- ymd("2023-10-31")

# Split data
btc_train <- btc_full %>% filter(`Open Time` <= split_date)
btc_test  <- btc_full %>% filter(`Open Time` > split_date)
eth_train <- eth_full %>% filter(`Open Time` <= split_date)
eth_test  <- eth_full %>% filter(`Open Time` > split_date)

# Helper function for returns at different frequencies
aggregate_returns <- function(df, time_col, freq) {
  df <- df %>%
    mutate(Time = lubridate::floor_date(!!sym(time_col), unit = freq)) %>%
    group_by(Time) %>%
    summarise(Close = dplyr::last(Close[!is.na(Close)]), .groups = "drop") %>%
    drop_na()
  
  ret <- TTR::ROC(df$Close, type = "discrete")[-1]
  df_ret <- data.frame(Time = df$Time[-1], Return = ret)
  return(df_ret)
}


# Frequencies to compute
frequencies <- c("day", "6 hours", "hour")

# Compute for BTC and ETH, full / train / test
returns_model_fitting <- list(
  btc_full = lapply(frequencies, function(f) aggregate_returns(btc_full, "Open Time", f)),
  btc_train = lapply(frequencies, function(f) aggregate_returns(btc_train, "Open Time", f)),
  btc_test = lapply(frequencies, function(f) aggregate_returns(btc_test, "Open Time", f)),
  eth_full = lapply(frequencies, function(f) aggregate_returns(eth_full, "Open Time", f)),
  eth_train = lapply(frequencies, function(f) aggregate_returns(eth_train, "Open Time", f)),
  eth_test = lapply(frequencies, function(f) aggregate_returns(eth_test, "Open Time", f))
)

names(returns_model_fitting$btc_full) <- frequencies
names(returns_model_fitting$btc_train) <- frequencies
names(returns_model_fitting$btc_test) <- frequencies
names(returns_model_fitting$eth_full) <- frequencies
names(returns_model_fitting$eth_train) <- frequencies
names(returns_model_fitting$eth_test) <- frequencies


# Select asset, select frequency, select, subset (train, test, full sample)
ret <- returns_model_fitting$eth_train$day
ret <- as.xts(ret)


max_p <- 1
max_q <- 1
max_r <- 3
max_s <- 3

plot_list <- list()


results_table <- data.frame(
  Series = character(),
  LLH = numeric(),
  AIC = numeric(),
  BIC = numeric(),
  Model = character(),
  stringsAsFactors = FALSE
)
test_table <- data.frame(
  Series = character(),
  JB = numeric(),
  BP = numeric(),
  ML = numeric(),
  LB = numeric(),
  stringsAsFactors = FALSE
)
aic_matrices <- list()
bic_matrices <- list()
llh_matrices <- list()

best_models <- list()


for (col in colnames(ret)) {
  series <- ret[, col]
  
  
  best_model <- NULL
  best_aic <- Inf
  best_bic <- Inf
  best_ll <- -Inf
  best_order <- NULL
  aic_matrix <- matrix(NA, nrow = max_p + 1, ncol = max_q + 1)
  bic_matrix <- matrix(NA, nrow = max_p + 1, ncol = max_q + 1)
  llh_matrix <- matrix(NA, nrow = max_p + 1, ncol = max_q + 1)
  jb_val <- NA
  bp_val <- NA
  ml_val <- NA
  lb_val <- NA
  
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      for (r in 0:max_r) {
        for (s in 0:max_s) {
          
          if (r == 0 && s == 0) {
            next
          }
          
          
          spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(r, s)),
                             mean.model = list(armaOrder = c(p, q)),
                             distribution.model = "std")
          fit <- tryCatch(ugarchfit(spec, series), error = function(e) NULL)
          
          if (!is.null(fit)) {
            aic <- -2*fit@fit$LLH + length(fit@fit$coef)
            bic <- -2*fit@fit$LLH + length(fit@fit$coef)*log(length(series))
            ll <- fit@fit$LLH
            
            aic_matrix[p + 1, q + 1] <- aic
            bic_matrix[p + 1, q + 1] <- bic
            llh_matrix[p + 1, q + 1] <- ll
            
            order <- paste("ARMA-GARCH(", p, ",", q, ";", r, ",", s, ")", sep = "")
            
            if (aic < best_aic) {
              best_model <- fit
              best_aic <- aic
              best_bic <- bic
              best_ll <- ll
              best_order <- order
              jb_val <- jarque.bera.test(best_model@fit$z)$statistic
              bp_val <- Box.test(best_model@fit$z^2, lag = 20, type = "Box-Pierce")$p.value
              ml_val <- Box.test(best_model@fit$z^2, lag = 30, type = "Ljung-Box")$p.value
              lb_val <- Box.test(best_model@fit$z, lag = 30, type = "Ljung-Box")$p.value
            }
          }
        }
      }
    }
  }
  
  aic_matrices[[col]] <- aic_matrix
  bic_matrices[[col]] <- bic_matrix
  llh_matrices[[col]] <- llh_matrix
  
  best_models[[col]] <- best_model
  
  results_table <- rbind(results_table, data.frame(
    Series = col,
    LLH = best_ll,
    AIC = best_aic,
    BIC = best_bic,
    Model = best_order
  ))
  test_table <- rbind(test_table, data.frame(
    Series = col,
    JB = jb_val,
    BP = bp_val,
    ML = ml_val,
    LB = lb_val
  ))
}
results_table


best_fit <- best_models[["BTC"]]
if (!is.null(best_fit)) {
  coef(best_fit)
  plot(best_fit, which = 2)
  plot(best_fit, which = 3)
} else {
  print("No valid model was fit for BTC.")
}
# ------------------------------
# GARCH Fitting
# ------------------------------

# Select frequ and asset
ret <- returns_model_fitting$eth_train$`6 hours`
ret <- as.xts(ret)


# Simple example using an ARMA (1,0) GARCH (1,1) for the BTC return series
model <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(2, 2)),
                    mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = "std")
fit <- ugarchfit(spec = model, data = ret)

# Params
print(fit)

# Output table
bic <- -2*fit@fit$LLH + length(fit@fit$coef)*log(length(ret))
aic <- -2*fit@fit$LLH + length(fit@fit$coef)
llh <- fit@fit$LLH
ml_val <- Box.test(fit@fit$z^2, lag = 30, type = "Ljung-Box")$p.value
lb_val <- Box.test(fit@fit$z, lag = 30, type = "Ljung-Box")$p.value



# Output plots: 2x2 grid
resid_std <- fit@fit$z
df <- coef(fit)["shape"]

plot(fit@fit$z)
qqplot(qt(ppoints(length(resid_std)), df = df), resid_std,
       main = "Q-Q Plot vs Student-t Distribution",
       xlab = "Theoretical Quantiles (t-distribution)",
       ylab = "Standardized Residuals")
acf(fit@fit$z)
acf(fit@fit$z^2)

# Looping over several model types and frequencies

# Define model types
model_types <- c("sGARCH", "eGARCH", "fGARCH", "gjrGARCH", "MS-GARCH")
frequencies <- names(returns_model_fitting$btc_train)

# Initialize AIC results table
bic_results <- data.frame(matrix(NA, nrow = length(frequencies), ncol = length(model_types)))
rownames(bic_results) <- frequencies
colnames(bic_results) <- model_types

# Loop through each frequency
for (freq in frequencies) {
  cat("Processing", freq, "\n")
  ret <- returns_model_fitting$btc_train[[freq]]
  ret <- as.xts(ret)
  
  for (model_type in model_types[1:4]) {
    if (model_type == "fGARCH") {
      spec <- ugarchspec(
        variance.model = list(model = model_type, garchOrder = c(1, 1), submodel = "TGARCH"),
        mean.model = list(armaOrder = c(1, 0)),
        distribution.model = "std"
      )
    } else {
      spec <- ugarchspec(
        variance.model = list(model = model_type, garchOrder = c(1, 1)),
        mean.model = list(armaOrder = c(1, 0)),
        distribution.model = "std"
      )
    }
    
    fit <- tryCatch(ugarchfit(spec = spec, data = ret, solver = "hybrid"),
                    error = function(e) NULL)
    
    if (!is.null(fit)) {
      bic_val <- -2*fit@fit$LLH + length(fit@fit$coef)*log(length(ret))
      bic_results[freq, model_type] <- bic_val
    }
  }
  
  
  ## MSGARCH model
  # Fit AR(1) model to get residuals for conditional variance estimation
  ar_fit <- tryCatch(arima(ret, order = c(1, 0, 0), include.mean = TRUE),
                     error = function(e) NULL)
  if (!is.null(ar_fit)) {
    ar_resid <- residuals(ar_fit)
    
    ms_spec <- CreateSpec(
      variance.spec = list(model = "eGARCH"),
      distribution.spec = list(distribution = "std"),
      switch.spec = list(do.mix = FALSE, K = 2)
    )
    
    ms_fit <- tryCatch(FitML(spec = ms_spec, data = ar_resid),
                       error = function(e) NULL)
    
    if (!is.null(ms_fit)) {
      llh <- ms_fit$loglik
      k <- length(ms_fit$par)
      n <- length(ms_fit$data)
      
      bic_val <- -2 * llh + log(n) * k
      bic_results[freq, "MS-GARCH"] <- bic_val
    }
  }
}

# View AIC results
print(bic_results)
btc_BICresults <- bic_results


# ------------------------------
# Rolling 1-step-ahead GARCH Forecasts
# ------------------------------

ret_full <- returns_model_fitting$btc_full$day
ret <- xts(ret_full$Return, order.by = ret_full$Time)
split_date <- as.Date("2023-10-31")

# Filter index of first prediction day (t+1)
start_index <- which(index(ret) > split_date)[1]

# Initialize storage
forecast_dates <- index(ret)[start_index:length(ret)]
forecast_values <- numeric(length(forecast_dates))
arch_values <- numeric(length(forecast_dates))   # alpha1
garch_values <- numeric(length(forecast_dates))  # beta1

# Model specification
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)

# Rolling forecast
for (i in seq_along(forecast_dates)) {
  idx <- start_index + i - 1
  ret_window <- ret[1:(idx - 1)]
  
  fit <- tryCatch(
    ugarchfit(spec = garch_spec, data = ret_window, solver = "hybrid"),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    forecast <- ugarchforecast(fit, n.ahead = 1)
    forecast_values[i] <- sigma(forecast)[1]
    
    coef_vec <- coef(fit)
    arch_values[i] <- as.numeric(fit@fit$coef["alpha1"])
    garch_values[i] <- as.numeric(fit@fit$coef["beta1"])
  } else {
    forecast_values[i] <- NA
    arch_values[i] <- NA
    garch_values[i] <- NA
  }
}

# Combine all into a single dataframe
garch_forecast_df <- data.frame(
  Date = forecast_dates,
  Forecast = forecast_values,
  ARCH = arch_values,
  GARCH = garch_values
)

ggplot(garch_forecast_df, aes(x = Date, y = ARCH)) +
  geom_line(color = "blue") +
  labs(
    title = "Rolling ARCH Parameter (α₁) Over Time",
    x = "Date",
    y = "ARCH Parameter (α₁)"
  ) +
  theme_minimal()

ggplot(garch_forecast_df, aes(x = Date, y = GARCH)) +
  geom_line(color = "red") +
  labs(
    title = "Rolling GARCH Parameter (β₁) Over Time",
    x = "Date",
    y = "GARCH Parameter (β₁)"
  ) +
  theme_minimal()




print(head(garch_forecast_df))
write.csv(garch_forecast_df, "forecasts/rolling_forecasts_btc_daily.csv", row.names = FALSE)

# Using the build in function

forecast.length <- sum(index(ret) > split_date)
window.size <- sum(index(ret) <= split_date)

# Run rolling forecast
roll <- ugarchroll(
  spec = garch_spec,
  data = ret,
  n.ahead = 1,
  forecast.length = forecast.length,
  refit.every = 7,
  refit.window = "moving",
  window.size = window.size,
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05),
  solver = "hybrid"
)

roll_df <- as.data.frame(roll)
print(head(roll_df))
write.csv(roll_df, "forecasts/rolling_forecasts_btc_daily.csv", row.names = FALSE)


# ------------------------------
# MS GARCH Fitting
# ------------------------------
ret <- returns_model_fitting$btc_train$day
ret <- as.xts(ret)

# Specify a two-regime sGARCH model with Student-t errors
spec <- CreateSpec(
  variance.spec = list(model = "eGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec = list(do.mix = FALSE, K = 2)  # 2 regimes, not a mixture model
)

# Fitting conditional mean
ar_fit <- arima(ret, order = c(1, 0, 0), include.mean = TRUE)
ar_resid <- residuals(ar_fit)

# Fitting conditional variance model
fit <- FitML(spec = spec, data = ar_resid)
fc <- predict(fit, nahead = 100, do.return.draw = FALSE)
fit$par
# Params and matrix:
fit

# Extract log-likelihood
llh <- fit$loglik

# Extract number of parameters
k <- length(fit$par)

# Sample size
n <- length(fit$data)

# Compute AIC and BIC
aic <- -2 * llh + 2 * k
bic <- -2 * llh + log(n) * k

# Extract standardized residuals

sigma_t_ms <- Volatility(fit)
z_ms <- as.numeric(ar_resid) / as.numeric(sigma_t_ms)

ml_val <- Box.test(z_ms^2, lag = 30, type = "Ljung-Box")$p.value
lb_val <- Box.test(z_ms, lag = 30, type = "Ljung-Box")$p.value


plot(z_ms)
acf(z_ms)
acf(z_ms^2)


# ------------------------------
# Rolling 1-step-ahead MS-eGARCH Forecasts (residuals from AR(1))
# ------------------------------
ret <- returns_model_fitting$btc_full$day
ret <- as.xts(ret)

split_date <- as.Date("2023-10-31")

# Filter index of first prediction day (t+1)
start_index <- which(index(ret) > split_date)[1]

# Initialize forecast storage
forecast_dates <- index(ret)[start_index:length(ret)]
forecast_values <- numeric(length(forecast_dates))

# Define MS-eGARCH(1,1) specification with 2 regimes
spec <- CreateSpec(
  variance.spec = list(model = "eGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec = list(K = 2, do.mix = FALSE)
)

# Initialize storage
sigma <- numeric(length(forecast_dates))

# Rolling forecast loop
for (i in seq_along(forecast_dates)) {
  idx <- start_index + i - 1
  ret_window_xts <- ret[1:(idx - 1)]  # subset xts
  ret_window <- as.numeric(ret_window_xts)
  
  # Step 1: Fit AR(1) with intercept to mean
  ar_fit <- tryCatch(arima(ret_window, order = c(1, 0, 0), include.mean = TRUE), error = function(e) NULL)
  
  if (!is.null(ar_fit)) {
    ar_resid <- residuals(ar_fit)
    
    # Step 2: Fit MS-eGARCH to residuals
    fit <- tryCatch(FitML(spec = spec, data = ar_resid), error = function(e) NULL)
    
    if (!is.null(fit)) {
      fc <- predict(fit, nahead = 1, do.return.draw = FALSE)
      sigma[i] <- as.numeric(fc$vol[1])
    } else {
      sigma[i] <- NA
    }
  } else {
    sigma[i] <- NA
  }
}

# Create output
msGarch_forecasts <- data.frame(
  Date = forecast_dates,
  Sigma = sigma
)


print(head(msGarch_forecasts))
write.csv(msGarch_forecasts, "forecasts/msgarch_forecasts_btc_daily.csv", row.names = FALSE)

# ------------------------------
# Filtered Historical Simulation: AR-eGARCH (Univariate)
# ------------------------------

# Prepare returns
ret <- returns_model_fitting$btc_full$day
ret_xts <- xts(ret$Return, order.by = ret$Time)

# Fit AR(1)-eGARCH(1,1) with Student-t errors
spec_uni <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)
fit_uni <- ugarchfit(spec = spec_uni, data = ret_xts)

# Extract standardized residuals
z_uni <- residuals/sigma(fit_uni)

# Forecast conditional mean and volatility
fc_uni <- ugarchforecast(fit_uni, n.ahead = 1)
mu_uni <- fitted(fc_uni)[1]
sigma_uni <- sigma(fc_uni)[1]

# Simulate future returns
B <- 10000
z_star_uni <- sample(z_uni, B, replace = TRUE)
r_sim_uni <- mu_uni + sigma_uni * z_star_uni

# Compute risk measures
VaR_99_uni <- quantile(r_sim_uni, 0.01)
VaR_95_uni <- quantile(r_sim_uni, 0.05)

cat("AR-eGARCH VaR 1%:", round(VaR_99_uni, 4), "\n")
cat("AR-eGARCH VaR 5%:", round(VaR_95_uni, 4), "\n")

# Plot simulated return distribution
hist(r_sim_uni, breaks = 50, col = "skyblue", main = "FHS: AR-eGARCH", xlab = "Simulated Return")
abline(v = VaR_99_uni, col = "red", lty = 2)
abline(v = VaR_95_uni, col = "orange", lty = 2)

# ------------------------------
# Filtered Historical Simulation: AR(1) + MS-eGARCH
# ------------------------------

# Prepare returns
ret <- returns_model_fitting$btc_full$day
ret_xts <- xts(ret$Return, order.by = ret$Time)
ret_num <- as.numeric(ret_xts)

# Fit AR(1) to returns
ar_fit <- arima(ret_num, order = c(1, 0, 0), include.mean = TRUE)
ar_resid <- residuals(ar_fit)

# Fit MS-eGARCH(1,1) with Student-t errors on residuals
spec_ms <- CreateSpec(
  variance.spec = list(model = "eGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec = list(K = 2, do.mix = FALSE)
)
fit_ms <- FitML(spec = spec_ms, data = ar_resid)

# Extract standardized residuals
sigma_t_ms <- Volatility(fit_ms)
z_ms <- ar_resid / sigma_t_ms

# Forecast conditional mean from AR(1)
phi <- ar_fit$coef["ar1"]
alpha <- if ("intercept" %in% names(ar_fit$coef)) ar_fit$coef["intercept"] else ar_fit$coef["mean"]
r_T <- tail(ret_num, 1)
mu_ms <- alpha + phi * r_T

# Forecast conditional volatility from MS-eGARCH
fc_ms <- predict(fit_ms, nahead = 1, do.return.draw = FALSE)
sigma_ms <- as.numeric(fc_ms$vol)

# Simulate future returns
B <- 10000
z_star_ms <- sample(z_ms, B, replace = TRUE)
r_sim_ms <- mu_ms + sigma_ms * z_star_ms

# Compute risk measures
VaR_99_ms <- quantile(r_sim_ms, 0.01)
VaR_95_ms <- quantile(r_sim_ms, 0.05)

cat("MS-eGARCH VaR 1%:", round(VaR_99_ms, 4), "\n")
cat("MS-eGARCH VaR 5%:", round(VaR_95_ms, 4), "\n")

# Plot simulated return distribution
hist(r_sim_ms, breaks = 50, col = "skyblue", main = "FHS: AR(1) + MS-eGARCH", xlab = "Simulated Return")
abline(v = VaR_99_ms, col = "red", lty = 2)
abline(v = VaR_95_ms, col = "orange", lty = 2)


# ------------------------------
# Historical Simulation: AR(1) + MS-eGARCH using simulate()
# ------------------------------

# Convert return series to numeric vector
ret_num <- as.numeric(ret_xts)

# MS eGARCH model with Student-t errors
spec <- CreateSpec(
  variance.spec = list(model = "eGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec = list(do.mix = FALSE, K = 2)
)

# Step 1: Fit AR(1) with intercept to mean
ar_fit <- arima(ret_num, order = c(1, 0, 0), include.mean = TRUE)
ar_resid <- residuals(ar_fit)

# Step 2: Fit MS-eGARCH to residuals
fit <- FitML(spec = spec, data = ar_resid)
fit$
sim <- simulate(object = fit, nsim = 1L, nahead = 10L,
                nburn = 500L)
head(sim)
plot(sim)
