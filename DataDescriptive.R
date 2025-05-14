library(xts)
library(tidyverse)
library(lubridate)
library(zoo)
library(tseries)
library(TTR)
library(forecast)
library(PerformanceAnalytics)
library(FinTS)
library(moments)
library(dplyr)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

# Load close price series
load_close_series <- function(path) {
  df <- read_csv(path, col_types = cols()) %>%
    select(`Open Time`, Close) %>%
    mutate(`Open Time` = ymd_hms(`Open Time`)) %>%
    arrange(`Open Time`) %>%
    drop_na()
  return(df)
}

# Compute descriptive statistics
compute_stats <- function(r, name) {
  vol_ann <- sd(r) * sqrt(252 * (1440 / name)) * 100
  return(c(
    "Average Return (%)" = mean(r) * 100,
    "Annualized Volatility (%)" = vol_ann,
    "Skewness" = skewness(r),
    "Excess Kurtosis" = kurtosis(r) - 3
  ))
}

# Plot time series
plot_series <- function(series_list, title, ylab) {
  df <- bind_rows(lapply(names(series_list), function(name) {
    data.frame(Date = index(series_list[[name]]), Value = coredata(series_list[[name]]), Asset = name)
  }))
  ggplot(df, aes(x = Date, y = Value, color = Asset)) +
    geom_line(alpha = 0.6) +
    labs(title = title, y = ylab, x = "Date") +
    theme_minimal()
}

# Daily annualized volatility from intraday returns
daily_ann_vol <- function(prices) {
  ret_1m <- ROC(prices$Close, type = "discrete")
  dates <- date(prices$`Open Time`[-1])
  df <- data.frame(Date = dates, Ret = ret_1m[-1]) %>%
    group_by(Date) %>%
    summarise(std = sd(Ret, na.rm = TRUE)) %>%
    mutate(vol_ann = std * sqrt(1440 * 252) * 100)
  return(df)
}

# Plot autocorrelogram
plot_autocorrelogram <- function(series, lags = 20, title = "", squared = FALSE) {
  if (squared) series <- series^2
  acf(series, lag.max = lags, main = title)
}

# Diagnostics tests
test_distribution_diagnostics <- function(series_list, lags = 20) {
  results <- lapply(names(series_list), function(name) {
    series <- na.omit(series_list[[name]])
    squared <- series^2
    lb <- Box.test(series, lag = lags, type = "Ljung-Box")
    ml <- Box.test(squared, lag = lags, type = "Ljung-Box")
    jb <- jarque.bera.test(series)
    data.frame(
      Series = name,
      Ljung_Box_Stat = lb$statistic,
      Ljung_Box_p = lb$p.value,
      McLeod_Li_Stat = ml$statistic,
      McLeod_Li_p = ml$p.value,
      Jarque_Bera_Stat = jb$statistic,
      Jarque_Bera_p = jb$p.value
    )
  })
  bind_rows(results)
}

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# ------------------------------
# Load and preprocess data
# ------------------------------

btc_full <- load_close_series("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/BTCUSDT_1m.csv")

eth_full <- load_close_series("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/ETHUSDT_1m.csv")

split_date <- ymd("2023-10-31")

# Split and aggregate prices at different frequencies
btc_train <- btc_full %>% filter(`Open Time` <= split_date)
eth_train <- eth_full %>% filter(`Open Time` <= split_date)

# Daily
btc_train_daily <- btc_train %>%
  mutate(Date = date(`Open Time`)) %>%
  group_by(Date) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()

eth_train_daily <- eth_train %>%
  mutate(Date = date(`Open Time`)) %>%
  group_by(Date) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()

# 1-hour
btc_train_hourly <- btc_train %>%
  mutate(Hour = floor_date(`Open Time`, unit = "hour")) %>%
  group_by(Hour) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()

eth_train_hourly <- eth_train %>%
  mutate(Hour = floor_date(`Open Time`, unit = "hour")) %>%
  group_by(Hour) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()

# 6-hour
btc_train_6hourly <- btc_train %>%
  mutate(SixHour = floor_date(`Open Time`, unit = "6 hours")) %>%
  group_by(SixHour) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()

eth_train_6hourly <- eth_train %>%
  mutate(SixHour = floor_date(`Open Time`, unit = "6 hours")) %>%
  group_by(SixHour) %>%
  summarise(Close = last(Close), .groups = "drop") %>%
  drop_na()


# ------------------------------
# Compute returns
# ------------------------------

# Daily
btc_ret_daily <- ROC(btc_train_daily$Close, type = "discrete")[-1]
eth_ret_daily <- ROC(eth_train_daily$Close, type = "discrete")[-1]

returns_daily_df <- data.frame(Date = btc_train_daily$Date[-1], BTC = btc_ret_daily, ETH = eth_ret_daily)
write_csv(returns_daily_df, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/train_returns_daily.csv")

# 1-hour
btc_ret_hourly <- ROC(btc_train_hourly$Close, type = "discrete")[-1]
eth_ret_hourly <- ROC(eth_train_hourly$Close, type = "discrete")[-1]

returns_hourly_df <- data.frame(Hour = btc_train_hourly$Hour[-1], BTC = btc_ret_hourly, ETH = eth_ret_hourly)
write_csv(returns_hourly_df, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/train_returns_hourly.csv")

# 6-hour
btc_ret_6hourly <- ROC(btc_train_6hourly$Close, type = "discrete")[-1]
eth_ret_6hourly <- ROC(eth_train_6hourly$Close, type = "discrete")[-1]

returns_6hourly_df <- data.frame(SixHour = btc_train_6hourly$SixHour[-1], BTC = btc_ret_6hourly, ETH = eth_ret_6hourly)
write_csv(returns_6hourly_df, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/train_returns_6hourly.csv")

# ------------------------------
# Descriptive statistics
# ------------------------------

# Daily
sum_stats_btc_ret_daily <- compute_stats(btc_ret_daily, 1440)
print(sum_stats_btc_ret_daily)

sum_stats_eth_ret_daily <- compute_stats(eth_ret_daily, 1440) 
print(sum_stats_eth_ret_daily)

# 1-hour
sum_stats_btc_ret_hourly <- compute_stats(btc_ret_hourly, 60)
print(sum_stats_btc_ret_hourly)

sum_stats_eth_ret_hourly <- compute_stats(eth_ret_hourly, 60) 
print(sum_stats_eth_ret_hourly)

# 6-hour
sum_stats_btc_ret_6hourly <- compute_stats(btc_ret_6hourly, 720)
print(sum_stats_btc_ret_6hourly)

sum_stats_eth_ret_6hourly <- compute_stats(eth_ret_6hourly, 720) 
print(sum_stats_eth_ret_6hourly)

# Print nice table with whole results
summary_table <- rbind(
  "BTC Daily" = sum_stats_btc_ret_daily,
  "BTC 6-Hour" = sum_stats_btc_ret_6hourly,
  "BTC Hourly" = sum_stats_btc_ret_hourly,
  "ETH Daily" = sum_stats_eth_ret_daily,
  "ETH 6-Hour" = sum_stats_eth_ret_6hourly,
  "ETH Hourly" = sum_stats_eth_ret_hourly
)

print(summary_table)

# ------------------------------
# Plots
# ------------------------------

# Prices

# Convert your data frames to XTS objects
# Daily
btc_daily_xts <- xts(btc_train_daily$Close, order.by = btc_train_daily$Date)
eth_daily_xts <- xts(eth_train_daily$Close, order.by = eth_train_daily$Date)

# 1-hour
btc_hourly_xts <- xts(btc_train_hourly$Close, order.by = btc_train_hourly$Hour)
eth_hourly_xts <- xts(eth_train_hourly$Close, order.by = eth_train_hourly$Hour)

# 6-hour
btc_6hourly_xts <- xts(btc_train_6hourly$Close, order.by = btc_train_6hourly$SixHour)
eth_6hourly_xts <- xts(eth_train_6hourly$Close, order.by = eth_train_6hourly$SixHour)

# Create named lists for the plot_series function
btc_daily_list <- list(BTC = btc_daily_xts)
eth_daily_list <- list(ETH = eth_daily_xts)
btc_hourly_list <- list(BTC = btc_hourly_xts)
eth_hourly_list <- list(ETH = eth_hourly_xts)
btc_6hourly_list <- list(BTC = btc_6hourly_xts)
eth_6hourly_list <- list(ETH = eth_6hourly_xts)

# Daily
plot_series(btc_daily_list, "Daily price series BTC", "BTC price USD")
plot_series(eth_daily_list, "Daily price series ETH", "ETH price USD")

# 1-hour
plot_series(btc_hourly_list, "Hourly price series BTC", "BTC price USD")
plot_series(eth_hourly_list, "Hourly price series ETH", "ETH price USD")

# 6-hour
plot_series(btc_6hourly_list, "6-Hour price series BTC", "BTC price USD")
plot_series(eth_6hourly_list, "6-Hour price series ETH", "ETH price USD")


# Returns

# Convert returns to xts objects
# Daily returns
btc_ret_daily_xts <- xts(btc_ret_daily, order.by = btc_train_daily$Date[-1])
eth_ret_daily_xts <- xts(eth_ret_daily, order.by = eth_train_daily$Date[-1])

# Hourly returns
btc_ret_hourly_xts <- xts(btc_ret_hourly, order.by = btc_train_hourly$Hour[-1])
eth_ret_hourly_xts <- xts(eth_ret_hourly, order.by = eth_train_hourly$Hour[-1])

# 6-hour returns
btc_ret_6hourly_xts <- xts(btc_ret_6hourly, order.by = btc_train_6hourly$SixHour[-1])
eth_ret_6hourly_xts <- xts(eth_ret_6hourly, order.by = eth_train_6hourly$SixHour[-1])

# Create named lists for the returns
btc_ret_daily_list <- list(BTC = btc_ret_daily_xts)
eth_ret_daily_list <- list(ETH = eth_ret_daily_xts)
btc_ret_hourly_list <- list(BTC = btc_ret_hourly_xts)
eth_ret_hourly_list <- list(ETH = eth_ret_hourly_xts)
btc_ret_6hourly_list <- list(BTC = btc_ret_6hourly_xts)
eth_ret_6hourly_list <- list(ETH = eth_ret_6hourly_xts)

# Daily
plot_series(btc_ret_daily_list, "Daily return series BTC", "BTC return")
plot_series(eth_ret_daily_list, "Daily return series ETH", "ETH return")

# 1-hour
plot_series(btc_ret_hourly_list, "Hourly return series BTC", "BTC return")
plot_series(eth_ret_hourly_list, "Hourly return series ETH", "ETH return")

# 6-hour
plot_series(btc_ret_6hourly_list, "6-Hour return series BTC", "BTC return")
plot_series(eth_ret_6hourly_list, "6-Hour return series ETH", "ETH return")


# ACF plots
# Daily
plot_autocorrelogram(btc_ret_daily_xts, lags = 20, title = "", squared = TRUE)
plot_autocorrelogram(eth_ret_daily_xts, lags = 20, title = "", squared = TRUE)

# 1-hour
plot_autocorrelogram(as.vector(btc_ret_hourly_xts), lags = 20, title = "", squared = TRUE)
plot_autocorrelogram(as.vector(eth_ret_hourly_xts), lags = 20, title = "", squared = TRUE)

# 6-hour
plot_autocorrelogram(as.vector(btc_ret_6hourly_xts), lags = 20, title = "", squared = TRUE)
plot_autocorrelogram(as.vector(eth_ret_6hourly_xts), lags = 20, title = "", squared = TRUE)


# PACF plots
# Daily


# 1-hour


# 6-hour


# Robust ACF plots
# Daily


# 1-hour


# 6-hour


# Robust PACF plots
# Daily


# 1-hour


# 6-hour

# ------------------------------
# Diagnostics
# ------------------------------





































