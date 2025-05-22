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
    mutate(Time = floor_date(!!sym(time_col), unit = freq)) %>%
    group_by(Time) %>%
    summarise(Close = last(Close), .groups = "drop") %>%
    drop_na()
  ret <- ROC(df$Close, type = "discrete")[-1]
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
ret <- returns_model_fitting$btc_train$day
ret <- as.xts(ret)


# Simple example using an ARMA (1,0) GARCH (1,1) for the BTC return series
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = "std")
fit <- ugarchfit(spec = model, data = ret)


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
            aic <- -2*fit@fit$LLH + 2*(r+s+p+q+1)
            bic <- -2*fit@fit$LLH + (r+s+p+q+1)*log(length(series))
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


# Select asset, select frequency, select, subset (train, test, full sample)
ret <- returns_model_fitting$btc_train$day
ret <- as.xts(ret)


# Univariate GARCH model forecasting
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = "std")
fit <- ugarchfit(spec = model, data = ret)




###
roll <- ugarchroll(
  spec, 
  data = returns,
  n.ahead = 1,
  forecast.length = 250,
  refit.every = 1,
  refit.window = "moving",
  window.size = 1000,
  solver = "hybrid",
  calculate.VaR = TRUE,
  VaR.alpha = c(0.01, 0.05)
)



# Univariate MS GARCH forecasting

# Specify a two-regime sGARCH model with Student-t errors
spec <- CreateSpec(
  variance.spec = list(model = "sGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec = list(do.mix = FALSE, K = 2)  # 2 regimes, not a mixture model
)

# Fit the model to your return series
fit <- FitML(spec = spec, data = ret)

# Show model summary
summary(fit)

# Plot smoothed regime probabilities
par(mfrow = c(2, 1))
plot(fit, which = "smoothed")   # Posterior regime probabilities
plot(fit, which = "conditional.sd")  # Regime-weighted conditional volatility

# One-step ahead volatility forecast (regime-weighted)
forecast <- Predict(spec, par = Coef(fit), data = ret, nahead = 1)
print(forecast)

plot(forecast)


sim <- Sim(spec = spec, par = Coef(fit), n = 1000)
plot(sim, type = "l", main = "Simulated MS-GARCH Returns")


# Return forecasting using FHS: Filtered Historical Simulation

# Univariate GARCH model

# Univariate MS GARCH model


# Filtered historical simulation
Fit ARMA-GARCH model and extract standardized residuals
Forecast next-period conditional mean and variance
Simulate future returns via bootstrapped standardized residuals
Compute quantiles / risk measures from simulated return distribution

