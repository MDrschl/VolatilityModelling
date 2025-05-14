library(tidyverse)
library(lubridate)
library(zoo)
library(tseries)
library(TTR)
library(forecast)
library(PerformanceAnalytics)
library(FinTS)
library(moments)

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

# Load and preprocess data
btc_full <- load_close_series("data/BTCUSDT_1m.csv")
eth_full <- load_close_series("data/ETHUSDT_1m.csv")
split_date <- ymd("2023-10-31")

# Split and aggregate daily prices
btc_train <- btc_full %>% filter(`Open Time` <= split_date)
eth_train <- eth_full %>% filter(`Open Time` <= split_date)

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

# Compute returns
btc_ret <- ROC(btc_train_daily$Close, type = "discrete")[-1]
eth_ret <- ROC(eth_train_daily$Close, type = "discrete")[-1]

returns_df <- data.frame(Date = btc_train_daily$Date[-1], BTC = btc_ret, ETH = eth_ret)
write_csv(returns_df, "data/train_returns.csv")
