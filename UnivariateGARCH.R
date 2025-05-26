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

par(mfrow=c(2,2))
acf(as.xts(returns_model_fitting$btc_train$day), main = "BTC Daily Returns")
acf(as.xts(returns_model_fitting$eth_train$day), main = "ETH Daily Returns")
acf(as.xts(returns_model_fitting$btc_train$day)^2, main = "BTC Daily Returns Squared")
acf(as.xts(returns_model_fitting$eth_train$day)^2, main = "ETH Daily Returns Squared")


compute_stats <- function(r, name) {
  vol_ann <- sd(r) * sqrt(252 * (1440 / name))
  return(c(
    "Mean" = round(mean(r),3),
    "Ann. Vola." = round(vol_ann,3),
    "Skew." = round(skewness(r), 3),
    "Excess K." = round(kurtosis(r) - 3,3), 
    "JB" = round(as.numeric(normalTest(r,method="jb")@test$statistic), 3)
  ))
}
# Daily
sum_stats_btc_ret_daily <- compute_stats(as.xts(returns_model_fitting$btc_train$day), 1440)
sum_stats_eth_ret_daily <- compute_stats(as.xts(returns_model_fitting$eth_train$day), 1440) 

# 6-hour
sum_stats_btc_ret_6hourly <- compute_stats(as.xts(returns_model_fitting$btc_train$`6 hours`), 720)
sum_stats_eth_ret_6hourly <- compute_stats(as.xts(returns_model_fitting$eth_train$`6 hours`), 720) 

# 1-hour
sum_stats_btc_ret_hourly <- compute_stats(as.xts(returns_model_fitting$btc_train$hour), 60)
sum_stats_eth_ret_hourly <- compute_stats(as.xts(returns_model_fitting$eth_train$hour), 60) 

# Print nice table with whole results
summary_table <- rbind(
  "BTC 24h" = sum_stats_btc_ret_daily,
  "BTC 6h" = sum_stats_btc_ret_6hourly,
  "BTC 1h" = sum_stats_btc_ret_hourly,
  "ETH 24h" = sum_stats_eth_ret_daily,
  "ETH 6h" = sum_stats_eth_ret_6hourly,
  "ETH 1h" = sum_stats_eth_ret_hourly
)

print(summary_table)
gamma=function(x,h)
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  gamma=sum(x[1:(n-h)]*x[(h+1):n])/n
}

rho=function(x,h)
{
  rho=gamma(x,h)/gamma(x,0)
}


gamma.asy=function(x,h) #Asymptotical cov. matrix of the sample autocorrelations
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  x2=x^2
  gamma.asy<-matrix(,h,h)  
  for (i in 1:h)
  {
    for (j in i:h)
    {
      gamma.asy[i,j]=gamma.asy[j,i]=sum(x[(j-i+1):(n-i)]*x[1:(n-j)]*x2[(j+1):n])/n
    }
  }
  rho.asy=1/gamma(x,0)^2*gamma.asy
  list(gamma.asy=gamma.asy,rho.asy=rho.asy)
}

#a=gamma.asy(microsoft,5)
corr.Box.test=function(x,h) #Corrected portmanteau test under GARCH assumption
{
  n<-length(x)
  a=gamma.asy(x,h)
  acf.val=sapply(c(1:h),function(h) rho(x,h))
  val=n*(acf.val%*%solve(a$rho.asy)%*%acf.val)
  p_val = 1-pchisq(val,h)
  return(c(val, p_val))
}
ret <- returns_model_fitting$btc_train$day
ret <- as.xts(ret)

Box.test(ret,lag=20,type="Ljung")
corr.Box.test(ret,20)


# Define metadata
assets <- c("btc", "eth")
frequencies <- c("day", "6 hours", "hour")
lags <- c(10, 20)
methods <- c("Q", "Q_star", "Q^2")  # Q = Ljung, Q_star = Corr, Q^2 = McLeod-Li

# Initialize results table
results_table <- expand.grid(
  Asset = assets,
  Frequency = frequencies,
  Lag = lags,
  stringsAsFactors = FALSE
)

# Initialize one column per method
for (m in methods) {
  results_table[[m]] <- NA_character_
}

# Updated test functions
corr.Box.test <- function(x, h) {
  n <- length(x)
  a=gamma.asy(x,h)
  acf.val <- sapply(c(1:h),function(h) rho(x,h))
  val <- n*(acf.val%*%solve(a$rho.asy)%*%acf.val)
  pval <- 1 - pchisq(val, h)
  return(c(stat = as.numeric(val), pval = as.numeric(pval)))
}

# Fill table
for (asset in assets) {
  for (freq in frequencies) {
    series_df <- returns_model_fitting[[paste0(asset, "_train")]][[freq]]
    series_xts <- xts(series_df$Return, order.by = series_df$Time)
    series_vec <- as.numeric(series_xts)
    
    for (lag in lags) {
      ljung_test <- Box.test(series_vec, lag = lag, type = "Ljung-Box")
      corr_test <- corr.Box.test(series_vec, lag)
      mcleod_test <- Box.test(series_vec^2, lag = lag, type = "Ljung-Box")
      
      # Format: statistic (p-value)
      Q_val <- sprintf("%.1f (%.2f)", ljung_test$statistic, ljung_test$p.value)
      Q_star_val <- sprintf("%.1f (%.2f)", corr_test["stat"], corr_test["pval"])
      Q2_val <- sprintf("%.1f (%.2f)", mcleod_test$statistic, mcleod_test$p.value)
      
      # Fill row
      row_idx <- which(
        results_table$Asset == asset &
          results_table$Frequency == freq &
          results_table$Lag == lag
      )
      results_table[row_idx, "Q"] <- Q_val
      results_table[row_idx, "Q_star"] <- Q_star_val
      results_table[row_idx, "Q^2"] <- Q2_val
    }
  }
}


results_wide <- pivot_wider(
  results_table,
  id_cols = c(Asset, Frequency),
  names_from = Lag,
  values_from = c(Q, Q_star, `Q^2`)
)

print(results_wide, row.names = FALSE)



# Select asset, select frequency, select, subset (train, test, full sample)
ret <- returns_model_fitting$btc_train
ret <- as.xts(ret)

max_p <- 1
max_q <- 1
max_r <- 3
max_s <- 3

results_table <- data.frame(
  Series = character(),
  BIC = numeric(),
  Model = character(),
  stringsAsFactors = FALSE
)


best_models <- list()


for (col in colnames(ret)) {
  series <- ret[, col]
  
  
  best_model <- NULL
  best_bic <- Inf
  best_order <- NULL

  
  for (r in 0:max_r) {
    for (s in 0:max_s) {
          
          if (r == 0 && s == 0) {
            next
          }
          
          
          spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(r, s)),
                             mean.model = list(armaOrder = c(1, 0)),
                             distribution.model = "std")
          fit <- tryCatch(ugarchfit(spec, series), error = function(e) NULL)
          
          if (!is.null(fit)) {
            llh <- fit@fit$LLH
            bic <- -2*llh + length(fit@fit$coef)*log(length(series))
            
            order <- paste("ARMA-GARCH(", 1, ",", 0, ";", r, ",", s, ")", sep = "")
            
            if (bic < best_bic) {
              best_model <- fit
              best_bic <- bic
              best_order <- order
            }
          }
      }
  }

  
  best_models[[col]] <- best_model
  
  results_table <- rbind(results_table, data.frame(
    Series = col,
    BIC = best_bic,
    Model = best_order
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
ret <- returns_model_fitting$btc_train$day
ret <- as.xts(ret)


# Simple example using an ARMA (1,0) GARCH (1,1) for the BTC return series
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = "std")
fit <- ugarchfit(spec = model, data = ret)
sqrt(diag(vcov(fit)))[4]
sqrt(diag(vcov(fit)))[5]

fit@fit$cvar
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

ret_full <- returns_model_fitting$btc_full$`6 hours`
ret <- xts(ret_full$Return, order.by = ret_full$Time)
#split_date <- as.Date("2023-10-31")
split_date <- as.Date("2025-01-01")

# Filter index of first prediction day (t+1)
start_index <- which(index(ret) > split_date)[1]

# Initialize storage
forecast_dates <- index(ret)[start_index:length(ret)]
forecast_values <- numeric(length(forecast_dates))

# Model specification
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)),
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
  } else {
    forecast_values[i] <- NA
  }
}
# Combine all into a single dataframe
garch_forecast_df <- data.frame(
  Date = forecast_dates,
  Forecast = forecast_values
)
print(head(garch_forecast_df))
write.csv(garch_forecast_df, "forecasts/sgarch_forecasts_btc_6h.csv", row.names = FALSE)

# Rolling forecast using only the last 800 observations
for (i in seq_along(forecast_dates)) {
  idx <- start_index + i - 1
  window_end <- idx - 1
  window_start <- max(1, window_end - 799)
  ret_window <- ret[window_start:window_end]
  
  fit <- tryCatch(
    ugarchfit(spec = garch_spec, data = ret_window, solver = "hybrid"),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    forecast <- ugarchforecast(fit, n.ahead = 1)
    forecast_values[i] <- sigma(forecast)[1]
  } else {
    forecast_values[i] <- NA
  }
}
# Combine all into a single dataframe
garch_forecast_df <- data.frame(
  Date = forecast_dates,
  Forecast = forecast_values
)

print(head(garch_forecast_df))
write.csv(garch_forecast_df, "forecasts/egarch_forecasts_btc_1h.csv", row.names = FALSE)


### Looping for param trajectories
ret_full <- returns_model_fitting$btc_full$day
ret <- xts(ret_full$Return, order.by = ret_full$Time)
split_date <- as.Date("2023-10-31")
#split_date <- as.Date("2025-01-01")

# Filter index of first prediction day (t+1)
start_index <- which(index(ret) > split_date)[1]

# Initialize storage
forecast_dates <- index(ret)[start_index:length(ret)]
forecast_values <- numeric(length(forecast_dates))

# Model specification
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)
# Initialize storage
forecast_values <- numeric(length(forecast_dates))
arch_values <- numeric(length(forecast_dates))
garch_values <- numeric(length(forecast_dates))
arch_lower <- numeric(length(forecast_dates))
arch_upper <- numeric(length(forecast_dates))
garch_lower <- numeric(length(forecast_dates))
garch_upper <- numeric(length(forecast_dates))

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

    alpha1 <- coef_vec["alpha1"]
    beta1  <- coef_vec["beta1"]
    se_alpha1 <- sqrt(diag(vcov(fit)))[4]
    se_beta1  <- sqrt(diag(vcov(fit)))[4]
    
    arch_values[i] <- alpha1
    garch_values[i] <- beta1
    
    arch_lower[i] <- alpha1 - 1.96 * se_alpha1
    arch_upper[i] <- alpha1 + 1.96 * se_alpha1
    garch_lower[i] <- beta1 - 1.96 * se_beta1
    garch_upper[i] <- beta1 + 1.96 * se_beta1
  } else {
    forecast_values[i] <- NA
    arch_values[i] <- NA
    garch_values[i] <- NA
    arch_lower[i] <- NA
    arch_upper[i] <- NA
    garch_lower[i] <- NA
    garch_upper[i] <- NA
  }
}
garch_forecast_df <- data.frame(
  Date = forecast_dates,
  Forecast = forecast_values,
  ARCH = arch_values,
  ARCH_Lower = arch_lower,
  ARCH_Upper = arch_upper,
  GARCH = garch_values,
  GARCH_Lower = garch_lower,
  GARCH_Upper = garch_upper
)

ggplot(garch_forecast_df, aes(x = Date)) +
  geom_ribbon(aes(ymin = ARCH_Lower, ymax = ARCH_Upper), fill = "skyblue", alpha = 0.3) +
  geom_line(aes(y = ARCH), color = "blue") +
  labs(
    title = "Rolling Estimate of ARCH Parameter with 95% CI",
    x = "Date",
    y = expression(alpha[1])
  ) +
  theme_minimal()

ggplot(garch_forecast_df, aes(x = Date)) +
  geom_ribbon(aes(ymin = GARCH_Lower, ymax = GARCH_Upper), fill = "salmon", alpha = 0.3) +
  geom_line(aes(y = GARCH), color = "red") +
  labs(
    title = "Rolling Estimate of GARCH Parameter with 95% CI",
    x = "Date",
    y = expression(beta[1])
  ) +
  theme_minimal()


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
ret <- returns_model_fitting$eth_train$day
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
print(fit)
# Extract log-likelihood
llh <- fit$loglik

# Extract number of parameters
k <- length(fit$par)

# Sample size
n <- length(fit$data)

# Compute AIC and BIC
aic <- -2 * llh + 2 * k
bic <- -2 * llh + log(n) * k

# Residuals from AR(1) mean model
ar_resid <- residuals(arima(ret, order = c(1, 0, 0), include.mean = TRUE))

# Volatility estimates (regime-specific)
sigma_t <- as.numeric(Volatility(fit))  # vector of length T

# Most likely regime at each t (1 or 2)
states <- State(fit)$Viterbi

# Extract estimated degrees of freedom
nu_1 <- fit$par["nu_1"]
nu_2 <- fit$par["nu_2"]

# Compute rescaling factor per t based on regime
scaling_factors <- ifelse(states == 1,
                          sqrt(nu_1 / (nu_1 - 2)),
                          sqrt(nu_2 / (nu_2 - 2)))

# Final standardized residuals
z_ms <- (as.numeric(ar_resid) / sigma_t) * scaling_factors

ml_val <- Box.test(z_ms^2, lag = 30, type = "Ljung-Box")$p.value
print(ml_val)
lb_val <- Box.test(z_ms, lag = 30, type = "Ljung-Box")$p.value
print(lb_val)

par(mfrow=c(1,1))
qqplot(qt(ppoints(length(z_ms)), df = df), z_ms,
       main = "Q-Q Plot vs Student-t Distribution",
       xlab = "Theoretical Quantiles (t-distribution)",
       ylab = "Standardized Residuals")
abline(0,1)
plot(z_ms)
acf(z_ms)
acf(z_ms^2)


# ------------------------------
# Rolling 1-step-ahead MS-eGARCH Forecasts (residuals from AR(1))
# ------------------------------
ret <- returns_model_fitting$eth_full$`6 hours`
ret <- as.xts(ret)

split_date <- as.Date("2025-01-01")

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
msGarch_forecasts_6h <- data.frame(
  Date = forecast_dates,
  Sigma = sigma
)


print(head(msGarch_forecasts_6h))
write.csv(msGarch_forecasts_6h, "forecasts/msgarch_forecasts_eth_6h.csv", row.names = FALSE)

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

sim <- simulate(object = fit, nsim = 1L, nahead = 100L,
                nburn = 500L)
head(sim)
plot(sim)
