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
library(ggplot2)

# Audrino's packages
library(fBasics)  # Load the package fBasics
library(zoo)
library(fGarch)   # Load the package for ARCH/GARCH estimation
library(tseries)
library(car)      # Consistent standard errors
library(systemfit)
library(mvtnorm)
library(quadprog)
library(VGAM)
library(sandwich)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

#-----------
# Load close price series
#-----------
load_close_series <- function(path) {
  df <- read_csv(path, col_types = cols()) %>%
    select(`Open Time`, Close) %>%
    mutate(`Open Time` = ymd_hms(`Open Time`)) %>%
    arrange(`Open Time`) %>%
    drop_na()
  return(df)
}

#-----------
# Compute descriptive statistics
#-----------
compute_stats <- function(r, name) {
  vol_ann <- sd(r) * sqrt(252 * (1440 / name)) * 100
  return(c(
    "Average Return (%)" = mean(r) * 100,
    "Annualized Volatility (%)" = vol_ann,
    "Skewness" = skewness(r),
    "Excess Kurtosis" = kurtosis(r) - 3
  ))
}

#-----------
# Plot time series
#-----------
plot_series <- function(series_list, title, ylab) {
  df <- bind_rows(lapply(names(series_list), function(name) {
    data.frame(Date = index(series_list[[name]]), Value = coredata(series_list[[name]]), Asset = name)
  }))
  ggplot(df, aes(x = Date, y = Value, color = Asset)) +
    geom_line(alpha = 0.6) +
    labs(title = title, y = ylab, x = "Date") +
    theme_minimal()
}

#-----------
# Daily annualized volatility from intraday returns
#-----------
daily_ann_vol <- function(prices) {
  ret_1m <- ROC(prices$Close, type = "discrete")
  dates <- date(prices$`Open Time`[-1])
  df <- data.frame(Date = dates, Ret = ret_1m[-1]) %>%
    group_by(Date) %>%
    summarise(std = sd(Ret, na.rm = TRUE)) %>%
    mutate(vol_ann = std * sqrt(1440 * 252) * 100)
  return(df)
}

#-----------
# Plot autocorrelogram (ACF)
#-----------
plot_autocorrelogram <- function(series, lags = 20, title = "", squared = FALSE) {
  if (squared) series <- series^2
  acf(series, lag.max = lags, main = title)
}

#-----------
# Plot partial autocorrelogram (PACF)
#-----------
plot_partial_autocorrelogram <- function(series, lags = 20, title = "", squared = FALSE) {
  if (squared) series <- series^2
  pacf(series, lag.max = lags, main = title)
}

#-----------
# Diagnostics tests
#-----------
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

#-----------
# Robust ACF plot
#-----------
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

n1.acf=function(x,main=NULL,method="NP")
{
  n=length(x)
  nlag=as.integer(min(10*log10(n),n-1))
  acf.val=sapply(c(1:nlag),function(h) rho(x,h))
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(acf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(acf.val,1.96*band,1.96/sqrt(n))
  acf(x,xlab="Lag",ylab="Sample autocorrelations",ylim=c(minval,maxval),main=main)
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}


#-----------
# PACF plot
#-----------
n1.pacf <- function(x, max_lag = 20, main = "Partial Autocorrelation Function") {
  # Calculate PACF values
  pacf_result <- pacf(x, lag.max = max_lag, plot = FALSE)
  pacf_values <- pacf_result$acf
  
  # Calculate standard confidence bands
  n <- length(x)
  conf_band <- 1.96/sqrt(n)
  
  # Create plot
  plot(1:max_lag, pacf_values, type = "h", 
       xlab = "Lag", ylab = "Partial Autocorrelation",
       main = main)
  
  # Add reference line at zero
  abline(h = 0, lty = 3)
  
  # Add confidence bands
  abline(h = conf_band, lty = 2, col = "blue")
  abline(h = -conf_band, lty = 2, col = "blue")
  
  # Return PACF values invisibly
  invisible(pacf_values)
}

#-----------
# RV + volatility signature plot
#-----------
volatility_signature <- function(data, 
                                 asset_name = "BTC", 
                                 start_date = "2023-10-23", 
                                 n_days = 20, 
                                 min_interval = 1, 
                                 max_interval = 60, 
                                 plot = TRUE){
  
  # ---------------
  # Data cleaning
  # ---------------
  
  day_date_list <- list()
  day_date_list[[1]] <- as.Date(start_date)
  
  for (d in 2:n_days){ 
    day_date_list[[d]] <- day_date_list[[d-1]] + 1
  }
  
  data_day <- list()
  data_day_clean <- list()
  data_day_formatted <- list()
  
  for (d in 1:n_days){
    
    data_day[[d]] <- subset(data, as.Date(`Open Time`) == day_date_list[[d]])
    
#    data_day_clean[[d]] <- data_day[[d]] %>%
#      group_by(`Open Time`) %>%
#      summarise(Close = last(Close), .groups = "drop")
    
    data_day_formatted[[d]] <- data_day[[d]] %>%
      mutate(
        Date = as.numeric(format(`Open Time`, "%Y%m%d")),
        Hour = as.numeric(format(`Open Time`, "%H")),
        Minute = as.numeric(format(`Open Time`, "%M")),
        Second = as.numeric(format(`Open Time`, "%S")),
        Volume = 1  # dummy value
      ) %>%
      select(Date, Hour, Minute, Second, Close, Volume)
  }
  
  # ---------------
  # RV calculation
  # ---------------
  
  sequence_list <- list()
  price_list <- list()
  log_returns_list <- list()
  rv_list <- list()
  
  for (d in 1:n_days){
    
    sequence_list[[d]] <- list() 
    price_list[[d]] <- list() 
    log_returns_list[[d]] <- list() 
    rv_list[[d]] <- numeric(max_interval - min_interval + 1)
    
    idx <- 1
    for (i in min_interval:max_interval){
      sequence_list[[d]][[idx]] <- seq(1, nrow(data_day_formatted[[d]]), by = i)      # List with vectors of indices
      price_list[[d]][[idx]] <- data_day_formatted[[d]][sequence_list[[d]][[idx]],]   # List with prices vectors 
      log_returns_list[[d]][[idx]] <- diff(log(price_list[[d]][[idx]][["Close"]]))    # List with log return vectors
      rv_list[[d]][idx] <- sum((log_returns_list[[d]][[idx]])^2)                      # List with estimate Realized Variance (at each interval)
      idx <- idx + 1
    }
  }
  
  # ---------------
  # Plot
  # ---------------
  
  if (plot) {
    # Set fixed y-axis limits to match the image
    ylim_values <- c(0.000, 0.001)
    
    plot(min_interval:max_interval, colMeans(do.call(rbind, rv_list)), type = "l", 
         main = paste("Volatility signature plot:", asset_name,
                      "\nFrom:", start_date, "for", n_days, "days"),
         xlab = "Minutes",
         ylab = "Sample RV",
         col = "black", 
         lwd = 1.5,
         ylim = ylim_values)  # Added ylim parameter here
    
    abline(h = mean(colMeans(do.call(rbind, rv_list))[floor(max_interval / 2):max_interval]), col = "grey", lty = "dotted")
  }
}

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# ------------------------------
# Load raw data inspection
# ------------------------------

# BTC
# Load raw data
btc_raw <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/BTCUSDT_1m.csv")

# Check structure and unique timestamps
str(btc_raw)
cat("Number of unique timestamps:", length(unique(btc_raw$Open.Time)), "\n")

# Count duplicates per timestamp (just to see)
dup_counts_btc <- table(btc_raw$Open.Time)
cat("Number of timestamps with duplicates:", sum(dup_counts_btc > 1), "\n")

# Clean data: aggregate duplicates (handles 2, 3, or more duplicates)
btc_raw_clean <- btc_raw %>%
  group_by(Open.Time) %>%
  summarise(
    Close = last(Close),
    .groups = "drop"
  ) %>%
  arrange(Open.Time)

cat("Cleaned data has", nrow(btc_raw_clean), "unique timestamp rows.\n")

# Optional: inspect first rows
print(head(btc_raw_clean[["Close"]]))

# Convert into format that we need
btc_full <- btc_raw_clean %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))

# ETH
# Load raw data
eth_raw <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/ETHUSDT_1m.csv")

# Check structure and unique timestamps
str(eth_raw)
cat("Number of unique timestamps:", length(unique(eth_raw$Open.Time)), "\n")

# Count duplicates per timestamp (just to see)
dup_counts_eth <- table(eth_raw$Open.Time)
cat("Number of timestamps with duplicates:", sum(dup_counts_eth > 1), "\n")

# Clean data: aggregate duplicates (handles 2, 3, or more duplicates)
eth_raw_clean <- eth_raw %>%
  group_by(Open.Time) %>%
  summarise(
    Close = last(Close),
    .groups = "drop"
  ) %>%
  arrange(Open.Time)

cat("Cleaned data has", nrow(eth_raw_clean), "unique timestamp rows.\n")

# Optional: inspect first rows
print(head(eth_raw_clean[["Close"]]))

# Convert into format that we need
eth_full <- eth_raw_clean %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))

# Check validity of minute data
summary(btc_full)
summary(eth_full)

# Save cleaned price series
write_csv(btc_full, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/btc_full.csv")

write_csv(eth_full, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/eth_full.csv")

# Plot minute level prices
# BTC plot


# ETH plot

# ------------------------------
# Preprocess data
# ------------------------------

# Load processed data
btc_full <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/btc_full.csv")
btc_full <- btc_full %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))

eth_full <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/eth_full.csv")
eth_full <- eth_full %>%
  rename(`Open Time` = `Open.Time`) %>%
  mutate(`Open Time` = ymd_hms(`Open Time`))

# Define split dayte
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

returns_hourly_df_2 <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/train_returns_hourly.csv")



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


# Robust ACF plots
# Daily
n1.acf(btc_ret_daily, main=c("BTC — daily"))
n1.acf(eth_ret_daily, main=c("ETH — daily"))


# 1-hour
n1.acf(btc_ret_hourly, main=c("BTC — hourly"))
n1.acf(eth_ret_hourly, main=c("ETH — hourly"))

# 6-hour
n1.acf(btc_ret_6hourly, main=c("BTC — 6 hourly"))
n1.acf(eth_ret_6hourly, main=c("ETH — 6 hourly"))


# PACF plots
# Daily
n1.pacf(btc_ret_daily, main=c("BTC — daily"))
n1.pacf(eth_ret_daily, main=c("ETH — daily"))

# 1-hour
n1.pacf(btc_ret_hourly, main=c("BTC — hourly"))
n1.pacf(eth_ret_hourly, main=c("ETH — hourly"))

# 6-hour
n1.pacf(btc_ret_6hourly, main=c("BTC — 6 hourly"))
n1.pacf(eth_ret_6hourly, main=c("ETH — 6 hourly"))

# ------------------------------
# Diagnostics
# ------------------------------

# Daily
diagnostics_btc_daily <- test_distribution_diagnostics(btc_daily_list, lags = 20)
diagnostics_eth_daily <- test_distribution_diagnostics(eth_daily_list, lags = 20)

# 1-hour
diagnostics_btc_hourly <- test_distribution_diagnostics(btc_hourly_list, lags = 20)
diagnostics_eth_hourly <- test_distribution_diagnostics(eth_hourly_list, lags = 20)

# 6-hour
diagnostics_btc_6hourly <- test_distribution_diagnostics(btc_6hourly_list, lags = 20)
diagnostics_eth_6hourly <- test_distribution_diagnostics(eth_6hourly_list, lags = 20)

# Print nice table with whole results
summary_table_diagnostics <- rbind(
  "BTC Daily" = diagnostics_btc_daily,
  "BTC 6-Hour" = diagnostics_btc_6hourly,
  "BTC Hourly" = diagnostics_btc_hourly,
  "ETH Daily" = diagnostics_eth_daily,
  "ETH 6-Hour" = diagnostics_eth_6hourly,
  "ETH Hourly" = diagnostics_eth_hourly
)

print(summary_table_diagnostics)

# ------------------------------
# Realized volatility
# ------------------------------

# Volatility signature plot

# BTC
# Start date = 18 August 2017, days = 14
volatility_signature(btc_full, 
                     asset_name = "BTC",
                     start_date = "2017-08-18",
                     n_days = 14,
                     min_interval = 1,
                     max_interval = 60,
                     plot = TRUE)

# Play around with starting date and days
# 2023
start_dates_2023 <- seq(as.Date("2023-01-02"), by = "7 days", length.out = 48)
n_days_seq = c(7, 14, 21, 28)


for (start_date in start_dates_2023){
  for (n_day in n_days_seq){
    volatility_signature(btc_full, 
                         asset_name = "BTC",
                         start_date = as.Date(start_date),
                         n_days = n_day,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# 2024
start_dates_2024 <- seq(as.Date("2024-01-02"), by = "7 days", length.out = 48)
n_days_seq = c(7, 14, 21, 28)

for (start_date in start_dates_2024){
  for (n_day in n_days_seq){
    volatility_signature(btc_full, 
                         asset_name = "BTC",
                         start_date = as.Date(start_date),
                         n_days = n_day,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# 2025
start_dates_2025 <- seq(as.Date("2025-01-02"), by = "7 days", length.out = 16)
n_days_seq = c(7, 14, 21, 28)

for (start_date in start_dates_2025){
  for (n_day in n_days_seq){
    volatility_signature(btc_full, 
                         asset_name = "BTC",
                         start_date = as.Date(start_date),
                         n_days = n_day,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# PLOTS USED
# Start date = 27 November 2023, days = 14
volatility_signature(btc_full, 
                     asset_name = "BTC",
                     start_date = "2023-11-27",
                     n_days = 14,
                     min_interval = 1,
                     max_interval = 60,
                     plot = TRUE)

volatility_signature(btc_full, 
                     asset_name = "BTC",
                     start_date = "2025-03-13",
                     n_days = 14,
                     min_interval = 1,
                     max_interval = 60,
                     plot = TRUE)



# ETH
# Play around with starting date and days
# 2017
start_dates_2017 <- seq(as.Date("2017-08-21"), by = "7 days", length.out = 1)
n_days_seq = c(7, 14, 21, 28)


for (start_date in start_dates_2017){
  for (n_day in n_days_seq){
    volatility_signature(eth_full, 
                         asset_name = "ETH",
                         start_date = as.Date(start_date),
                         n_days = 14,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# 2023
start_dates_2023 <- seq(as.Date("2023-01-02"), by = "7 days", length.out = 48)
n_days_seq = c(7, 14, 21, 28)


for (start_date in start_dates_2023){
  for (n_day in n_days_seq){
    volatility_signature(eth_full, 
                         asset_name = "ETH",
                         start_date = as.Date(start_date),
                         n_days = n_day,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# 2024
start_dates_2024 <- seq(as.Date("2024-01-02"), by = "7 days", length.out = 48)
n_days_seq = c(7, 14, 21, 28)


for (start_date in start_dates_2024){
  for (n_day in n_days_seq){
    volatility_signature(eth_full, 
                         asset_name = "ETH",
                         start_date = as.Date(start_date),
                         n_days = 14,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# 2025
start_dates_2025 <- seq(as.Date("2025-01-02"), by = "7 days", length.out = 16)
n_days_seq = c(7, 14, 21, 28)

for (start_date in start_dates_2025){
  for (n_day in n_days_seq){
    volatility_signature(eth_full, 
                         asset_name = "ETH",
                         start_date = as.Date(start_date),
                         n_days = n_day,
                         min_interval = 1,
                         max_interval = 60,
                         plot = TRUE)
  }
}

# Plots used
volatility_signature(eth_full, 
                     asset_name = "ETH",
                     start_date = "2023-09-25",
                     n_days = 14,
                     min_interval = 1,
                     max_interval = 60,
                     plot = TRUE)

volatility_signature(eth_full, 
                     asset_name = "ETH",
                     start_date = "2025-03-13",
                     n_days = 14,
                     min_interval = 1,
                     max_interval = 60,
                     plot = TRUE)


# Extract Realized Variances for out-of-sample testing
# BTC - 1 minute intervals
btc_rv <- volatility_signature(btc_full,
                               asset_name = "BTC",
                               start_date = "2023-11-01",
                               n_days = 546,
                               min_interval = 1,
                               max_interval = 1,
                               plot = FALSE)

rv_values_btc <- unlist(btc_rv$BTC)
rv_dates_btc <- seq.Date(from = as.Date("2023-11-01"), by = "day", length.out = 546)
btc_rv_df <- data.frame(Date = rv_dates_btc, RV_1min = rv_values_btc)

write_csv(btc_rv_df, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/btc_rv_df.csv")

# ETH - 1 minute intervals
eth_rv <- volatility_signature(eth_full, 
                               asset_name = "ETH",
                               start_date = "2023-11-01",
                               n_days = 547,
                               min_interval = 1,
                               max_interval = 1,
                               plot = FALSE)

rv_values_eth <- unlist(eth_rv$ETH)
rv_dates_eth <- seq.Date(from = as.Date("2023-11-01"), by = "day", length.out = 547)
eth_rv_df <- data.frame(Date = rv_dates_eth, RV_1min = rv_values_eth)

write_csv(eth_rv_df, "/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/eth_rv_df.csv")









































