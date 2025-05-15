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

# Audrino's packages
library(fBasics) #Load the package fBasics
library(zoo)
library(fGarch) #Load the package for ARCH/GARCH estimation
library(tseries)
library(car) #consistent standard errors
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
  nlag=40
  acf.val=sapply(c(1:nlag),function(h) rho(x,h))
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(acf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(acf.val,1.96*band,1.96/sqrt(n))
  acf(x,lag.max=nlag,xlab="Lag",ylab="Sample autocorrelations",ylim=c(minval,maxval),main=main)
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}

#-----------
# Work-around
#-----------
n2.acf = function(x, main = NULL, method = "NP")
{
  n = length(x)
  nlag = 40
  
  # Calculate autocorrelations manually 
  acf.val = sapply(c(1:nlag), function(h) rho(x, h))
  
  # Calculate robust bands
  x2 = x^2
  var = 1 + (sapply(c(1:nlag), function(h) gamma(x2, h))) / gamma(x, 0)^2
  band = sqrt(var/n)
  
  # Calculate regular confidence band for comparison
  std.band = 1.96/sqrt(n)
  
  # Set plot limits
  minval = 1.2 * min(acf.val, -1.96 * band, -std.band)
  maxval = 1.2 * max(acf.val, 1.96 * band, std.band)
  
  # Calculate ACF using R's function but don't plot it
  acf_obj = acf(x, lag.max = nlag, plot = FALSE)
  
  # Manual plotting
  plot(0:nlag, acf_obj$acf, type = "h", 
       xlab = "Lag", ylab = "Sample autocorrelations",
       ylim = c(minval, maxval), main = main)
  
  # Add reference line at zero
  abline(h = 0)
  
  # Add the robust confidence bands
  lines(1:nlag, 1.96 * band, lty = 1, col = "red")
  lines(1:nlag, -1.96 * band, lty = 1, col = "red")
  
  # Optionally add standard bands for comparison
  abline(h = std.band, lty = 2, col = "blue")
  abline(h = -std.band, lty = 2, col = "blue")
}

#-----------
# Simplified robust ACF function following Francq & Zakoian (2009)
#-----------
simple_robust_acf = function(x, main = NULL, lags = 30) {
  # Center the series
  x = x - mean(x)
  n = length(x)
  
  # Get the ACF values using acf() for convenience
  acf_obj = acf(x, lag.max = lags, plot = FALSE)
  acf_values = as.vector(acf_obj$acf[-1])  # Remove lag 0
  
  # Standard error for standard bands
  std_se = 1/sqrt(n)
  
  # Create the plot without confidence bands
  plot(0:lags, c(1, acf_values), type = "h",
       ylim = c(-0.2, 1),
       xlab = "Lag", ylab = "ACF", main = main)
  abline(h = 0)
  
  # Add standard confidence bands
  abline(h = 1.96 * std_se, lty = 2, col = "blue")
  abline(h = -1.96 * std_se, lty = 2, col = "blue")
  
  # Add simple scaled bands for GARCH effects (approximately 2x the standard)
  abline(h = 3.92 * std_se, col = "red")  # 2x the standard
  abline(h = -3.92 * std_se, col = "red") # 2x the standard
  
  cat("Standard 95% bounds: ±", round(1.96 * std_se, 4), "\n")
  cat("Robust 95% bounds: ±", round(3.92 * std_se, 4), "\n")
}

#-----------
# Robust PACF plot
#-----------
n1.pacf=function(x,main=NULL,method="NP")
{
  n=length(x)
  nlag=40
  pacf_obj = pacf(x, lag.max = nlag, plot = FALSE)
  pacf.val = pacf_obj$acf
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(pacf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(pacf.val,1.96*band,1.96/sqrt(n))
  pacf(x,lag.max=nlag,xlab="Lag",ylab="Sample partial autocorrelations",ylim=c(minval,maxval),main=main)
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}

#-----------
# Work-around PACF
#-----------
n2.pacf = function(x, main = NULL, method = "NP")
{
  n = length(x)
  nlag = 40
  
  # Calculate robust bands
  x2 = x^2
  var = 1 + (sapply(c(1:nlag), function(h) gamma(x2, h))) / gamma(x, 0)^2
  band = sqrt(var/n)
  
  # Calculate regular confidence band for comparison
  std.band = 1.96/sqrt(n)
  
  # Calculate PACF using R's function but don't plot it
  pacf_obj = pacf(x, lag.max = nlag, plot = FALSE)
  pacf.val = pacf_obj$acf
  
  # Set plot limits
  minval = 1.2 * min(pacf.val, -1.96 * band, -std.band)
  maxval = 1.2 * max(pacf.val, 1.96 * band, std.band)
  
  # Manual plotting
  plot(1:nlag, pacf.val, type = "h", 
       xlab = "Lag", ylab = "Sample partial autocorrelations",
       ylim = c(minval, maxval), main = main)
  
  # Add reference line at zero
  abline(h = 0)
  
  # Add the robust confidence bands
  lines(1:nlag, 1.96 * band, lty = 1, col = "red")
  lines(1:nlag, -1.96 * band, lty = 1, col = "red")
  
  # Add standard bands for comparison
  abline(h = std.band, lty = 2, col = "blue")
  abline(h = -std.band, lty = 2, col = "blue")
}

#-----------
# Function to compute the RV at different frequencies
#-----------
hfrtn <- function(da,int,logrtn=TRUE){
  # Compute intraday returns
  #
  # int: time intervals in minutes
  # da: data in the format: date, hour, minute, second, price, volume
  #
  if(!is.matrix(da))da=as.matrix(da)
  intsec=int*60
  istart=0
  iend=24*60*60  
  # compute the number of prices
  tradetime=iend - istart
  ntrade=floor(tradetime/intsec)
  T=dim(da)[1]
  nday=da[T,1]-da[1,1]+1
  npri=nday*ntrade
  #print(c(ntrade,nday,npri))
  
  price=rep(0,npri)
  # price is the last transaction price of the time interval
  caltime = da[,2]*60*60 + da[,3]*60
  #caltime=da[,2]*60*60+da[,3]*60+da[,4]
  #plot(caltime,type='l')
  
  icnt=0
  date=da[1,1]
  for (i in 1:T) {
    if(caltime[i] >= istart){
      iday=da[i,1]-date
      if(caltime[i] < (iend+1)){
        
        if(caltime[i]==iend){
          price[iday*ntrade+ntrade]=da[i,5]
        }
        
        if((caltime[i] >= istart) && (caltime[i] < iend)){
          ii=caltime[i]-istart
          ij=floor(ii/intsec)
          price[iday*ntrade+ij+1]=da[i,5]
        }
      }
    }
  }
  for (i in 2:npri){
    if(price[i] <= 0)price[i]=price[i-1]
  }
  
  plot(price,type='l')
  
  pri=log(price)
  #skip overnight returns
  nrtn=ntrade-1
  rtn=NULL
  for (i in 1:nday){
    ist=(i-1)*ntrade
    for (j in 2:ntrade){
      rtn=c(rtn,pri[ist+j]-pri[ist+j-1])
    }
  }
  
  hfrtn = list(rtn=rtn,price=price)
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
plot_partial_autocorrelogram(btc_ret_daily_xts, lags = 20, title = "", squared = FALSE)
plot_partial_autocorrelogram(eth_ret_daily_xts, lags = 20, title = "", squared = FALSE)

# 1-hour
plot_partial_autocorrelogram(as.vector(btc_ret_hourly_xts), lags = 20, title = "", squared = FALSE)
plot_partial_autocorrelogram(as.vector(eth_ret_hourly_xts), lags = 20, title = "", squared = FALSE)

# 6-hour
plot_partial_autocorrelogram(as.vector(btc_ret_6hourly_xts), lags = 20, title = "", squared = FALSE)
plot_partial_autocorrelogram(as.vector(eth_ret_6hourly_xts), lags = 20, title = "", squared = FALSE)


# Robust ACF plots
# Daily
n1.acf(btc_ret_daily_xts, main = c(""))
n1.acf(eth_ret_daily_xts, main = c(""))

# 1-hour
n2.acf(btc_ret_hourly_xts, main = c(""))
n2.acf(eth_ret_hourly_xts, main = c(""))

# 6-hour
n2.acf(btc_ret_6hourly_xts, main = c(""))
n2.acf(eth_ret_6hourly_xts, main = c(""))


# Robust PACF plots
# Daily
n1.pacf(btc_ret_daily_xts, main = c(""))
n1.pacf(eth_ret_daily_xts, main = c(""))

# 1-hour
n2.pacf(btc_ret_hourly_xts, main = c(""))
n2.pacf(eth_ret_hourly_xts, main = c(""))

# 6-hour
n2.pacf(btc_ret_6hourly_xts, main = c(""))
n2.pacf(eth_ret_6hourly_xts, main = c(""))

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

# Load the data in the correct format

# BTC
day1_date_btc <- as.Date("2023-10-23")
day2_date_btc <- day1_date_btc + 1
day3_date_btc <- day2_date_btc + 1
day4_date_btc <- day3_date_btc + 1
day5_date_btc <- day4_date_btc + 1

days_list <- list(day1_date_btc, day2_date_btc, day3_date_btc, day4_date_btc ,day5_date_btc)

btc_day <- list()
btc_day_clean <- list()
btc_day_formatted <- list()

for (d in 1:5){
  
  btc_day[[d]] <- subset(btc_full, as.Date(`Open Time`) == days_list[[d]])
  
  btc_day_clean[[d]] <- btc_day[[d]] %>%
    group_by(`Open Time`) %>%
    summarise(Close = last(Close), .groups = "drop")
  
  btc_day_formatted[[d]] <- btc_day_clean[[d]] %>%
    mutate(
      Date = as.numeric(format(`Open Time`, "%Y%m%d")),
      Hour = as.numeric(format(`Open Time`, "%H")),
      Minute = as.numeric(format(`Open Time`, "%M")),
      Second = as.numeric(format(`Open Time`, "%S")),
      Volume = 1  # dummy value
    ) %>%
    select(Date, Hour, Minute, Second, Close, Volume)
}

# Calculate realized variances
rv.btc=matrix(,5,60)

for (d in 1:5){
  for (t in 1:60)
  {
    a=hfrtn(da=btc_day_formatted[[d]],int=t)
    
    for (i in 1:length(a$rtn))
    {
      if (is.nan(a$rtn[i]))
      {
        a$rtn[i]=0
      }
    }
    
    rv.btc[d,t]=sum(a$rtn^2)
  }
}

rv=rv.btc
#a=rv
a=colMeans(rv)

ts.plot(a,xlab="Minutes",ylab="Sample RV",main="Volatility signature plot") #Volatility signature plot
lines(rep(mean(a[30:60]),60),lty=3)

# ETH


























