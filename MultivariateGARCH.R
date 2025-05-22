library(rugarch)
library(rmgarch)
library(zoo)
library(xts)

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


# ------------------------------
# Model selection based on lowest BIC
# ------------------------------


retBTC <- returns_model_fitting$btc_train$day
retETH <- returns_model_fitting$eth_train$day

retBTC_xts <- xts(retBTC$Return, order.by = as.Date(retBTC$Time))
retETH_xts <- xts(retETH$Return, order.by = as.Date(retETH$Time))
multiReturns <- merge(retBTC_xts, retETH_xts)
colnames(multiReturns) <- c("BTC", "ETH")

dcc_orders <- list(c(1, 0), c(0, 1), c(1, 1), c(2, 1), c(1, 2), c(2,2))
distributions <- c("mvnorm", "mvt")

results <- list()
best_bic <- Inf
best_model <- NULL
best_spec_desc <- ""

for (dcc_order in dcc_orders) {
  for (dist in distributions) {
    spec <- dccspec(
      uspec = multispec(list(
        ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(1, 0)),
                   distribution.model = "std"),
        ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(1, 0)),
                   distribution.model = "std"))),
      dccOrder = dcc_order,
      distribution = dist
    )
    fit <- tryCatch(dccfit(spec, data = multiReturns), error = function(e) NULL)
    
    if (!is.null(fit)) {
      bic_val <- infocriteria(fit)[2]
      spec_desc <- paste0("DCC(", dcc_order[1], ",", dcc_order[2], "), Dist: ", dist)
      results[[spec_desc]] <- list(fit = fit, bic = bic_val)
      
      if (bic_val < best_bic) {
        best_bic <- bic_val
        best_model <- fit
        best_spec_desc <- spec_desc
      }
    } else {
      spec_desc <- paste0("DCC(", dcc_order[1], ",", dcc_order[2], "), Dist: ", dist)
      results[[spec_desc]] <- "Fitting failed"
    }
  }
}

print("Results for all specifications:")
print(lapply(results, function(res) if (is.list(res)) paste("BIC:", round(res$bic, 3)) else res))

print(paste("\nBest model based on BIC:", best_spec_desc, "with BIC:", round(best_bic, 3)))
coef(best_model)


# ------------------------------
# Fitting best model
# ------------------------------


spec <- dccspec(
  uspec = multispec(list(
    ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
               mean.model = list(armaOrder = c(1, 0)),
               distribution.model = "std"),
    ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
               mean.model = list(armaOrder = c(1, 0)),
               distribution.model = "std"))),
  dccOrder = c(1,1),
  distribution = "mvt")

fit <- dccfit(spec, data = multiReturns)
print(fit)
plot(fit)
correlation <- rcor(fit)


# Extract dynamic correlation between BTC and ETH (off-diagonal)
btc_eth_corr <- correlation[1, 2, ]
btc_eth_corr <- xts(btc_eth_corr, order.by = index(multiReturns))

autoplot(as.zoo(btc_eth_corr)) +
  labs(title = "Dynamic Correlation between BTC and ETH", y = "Correlation") +
  theme_minimal()

grangertest(ETH^2 ~ BTC^2, order = 1, data = as.data.frame(multiReturns))

# Checking for vola spillover, but there might be more sophisticated checks for this
BTC_vol <- sigma(fit)[, 1]
ETH_vol <- sigma(fit)[, 2]

summary(lm(ETH_vol ~ lag(BTC_vol, 1)))

# ------------------------------
# Rolling 1-step-ahead DCC GARCH Forecasts (BTC & ETH)
# ------------------------------

retBTC <- returns_model_fitting$btc_full$day
retETH <- returns_model_fitting$eth_full$day

# Prepare data
retBTC_xts <- xts(retBTC$Return, order.by = as.Date(retBTC$Time))
retETH_xts <- xts(retETH$Return, order.by = as.Date(retETH$Time))
multiReturns <- merge(retBTC_xts, retETH_xts, join = "inner")
colnames(multiReturns) <- c("BTC", "ETH")

# Define forecast start
split_date <- as.Date("2023-10-31")
start_index <- which(index(multiReturns) > split_date)[1]

# Define forecast horizon and dates
forecast_dates <- index(multiReturns)[start_index:nrow(multiReturns)]
forecast_values_btc <- numeric(length(forecast_dates))
forecast_values_eth <- numeric(length(forecast_dates))

# DCC model specification
dcc_spec <- dccspec(
  uspec = multispec(list(
    ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    ),
    ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    )
  )),
  dccOrder = c(1, 1),
  distribution = "mvt"
)

# Rolling DCC-GARCH forecast loop
for (i in seq_along(forecast_dates)) {
  idx <- start_index + i - 1
  ret_window <- multiReturns[1:(idx - 1), ]
  
  fit <- tryCatch(
    dccfit(dcc_spec, data = ret_window, fit.control = list(eval.se = FALSE)),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    fc <- dccforecast(fit, n.ahead = 1)
    H_t1 <- rcov(fc)[[1]][, , 1]
    forecast_values_btc[i] <- sqrt(H_t1[1, 1])  # BTC volatility
    forecast_values_eth[i] <- sqrt(H_t1[2, 2])  # ETH volatility
  } else {
    forecast_values_btc[i] <- NA
    forecast_values_eth[i] <- NA
  }
}

# Compile result
dcc_forecast_df <- data.frame(
  Date = forecast_dates,
  BTC_StdDev = forecast_values_btc,
  ETH_StdDev = forecast_values_eth
)

# Output
print(head(dcc_forecast_df))
write.csv(dcc_forecast_df, "forecasts/rolling_dcc_forecast_btc_eth.csv", row.names = FALSE)



### Example for one step ahead

retBTC <- returns_model_fitting$btc_full$day
retETH <- returns_model_fitting$eth_full$day
retBTC_xts <- xts(retBTC$Return, order.by = as.Date(retBTC$Time))
retETH_xts <- xts(retETH$Return, order.by = as.Date(retETH$Time))
multiReturns <- merge(retBTC_xts, retETH_xts)
colnames(multiReturns) <- c("BTC", "ETH")

dcc_spec <- dccspec(
  uspec = multispec(list(
    ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
               mean.model = list(armaOrder = c(1, 0)),
               distribution.model = "std"),
    ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
               mean.model = list(armaOrder = c(1, 0)),
               distribution.model = "std")
  )),
  dccOrder = c(1, 1),
  distribution = "mvt"
)

fit = dccfit(dcc_spec, multiReturns[1:which(index(multiReturns) == split_date)])

split_date <- as.Date("2023-10-31")
forecasting_date <- as.Date("2023-11-5")

y = dccforecast(fit)


roll <- dccroll(
  spec = dcc_spec,
  data = multiReturns[1:which(index(multiReturns) == forecasting_date),],
  n.ahead = 1,
  n.start =  which(index(multiReturns) == split_date),
  refit.every = 1,
  refit.window = "moving",
  solver = "solnp",
  fit.control = list(eval.se = FALSE),
  solver.control = list(trace = 0)
)

# 5. Extract correlation matrix forecast for 2023-11-01
forecast_date <- index(multiReturns)[n_start + 1]
forecast_corr <- rcor(roll)[[1]][, , 1]  # 1-step ahead correlation

# Extract BTC-ETH correlation
btc_eth_corr <- forecast_corr[1, 2]

cat("Forecasted BTC–ETH correlation for", forecast_date, "is", round(btc_eth_corr, 4), "\n")
