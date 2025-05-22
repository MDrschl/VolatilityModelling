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
library(rugarch)


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

#-----------
# XXX
#-----------














# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# ------------------------------
# Align data
# ------------------------------

# Placeholder to be replaced by actual volatilits forecasts
vol_forecast <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_rv_df$realized_vol + rnorm(nrow(btc_rv_df), 0, 0.01)  # small noise
)

# Check length and specifications of data first
head(vol_forecast)
tail(vol_forecast)

head(btc_rv_df)
tail(btc_rv_df)

length(vol_forecast$Date) == length(btc_rv_df$Date)
all(vol_forecast$Date == btc_rv_df$Date)

# Convert realized variance sequence into realized volatility sequence
btc_rv_df$realized_vol <- sqrt(btc_rv_df$RV_1min)

# ------------------------------
# Visual comparison
# ------------------------------

# Scatter plot comparison
plot(vol_forecast$forecast_vol, btc_rv_df$realized_vol, xlab="XX GARCH(X,X) predictions",ylab="Volatility proxy")
abline(0,1,lty=3)

# Time-series plot overlap









# ------------------------------
# Visual comparison
# ------------------------------















