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
# SPA test
#-----------
spa <- function(per = perf, bench = 1, m = 4, obs = 546, q = 0.25, iter = 1000, periodogram = T){
  
  # Test of superior predictive ability of Hansen, 2005 JBES
  e = bench # benchmark
  d = matrix(, obs, m-1)
  s = 0
  
  for (i in seq(1, m, 1)[-e]){
    
    s = s+1
    d[, s] = per[, e] - per[, i]
    
  }
  
  #colMeans(d)
  
  w = rep(0, m-1)
  
  for (k in 1:(m-1)){
    
    #e=c()
    #for (i in 0:(obs-1))
    #{
    #e=c(e,1/obs*sum((d[1:(obs-i),k]-mean(d[,k]))*(d[(1+i):(obs),k]-mean(d[,k]))))
    e = acf(d[, k], lag.max = obs-1, type = "covariance", plot = F)$acf
    #}
    
    if (periodogram == F){
      w[k] = sqrt(e[1]+2 * sum(((obs-seq(1, obs-1, 1)) / obs*(1-q)^{seq(1, obs-1, 1)}+seq(1, obs-1, 1)/obs*(1-q)^{obs-seq(1, obs-1, 1)})*e[2:obs]))
    }
    
    else if (periodogram == T){
      w[k] = sqrt(spectrum(d[,k],plot = F)$spec[1])
    }
  }
  
  #print(sqrt(obs)*colMeans(d)/w)
  #print(pnorm(sqrt(obs)*colMeans(d)/w))
  
  stat = max(0,max(sqrt(obs)*colMeans(d)/w))
  
  #Bootstrap:
  
  stat.boos = rep(0,iter)
  
  for (r in 1:iter){
    
    #print(r)
    tau = rep(0,obs)
    tau[1] = as.integer(obs*runif(1))+1
    
    for (i in 2:obs){
      
      s = runif(1)
      tau[i] = (as.integer(obs*runif(1))+1)*(s<q)+((tau[i-1]<obs)*tau[i-1]+1)*(s >= q)
      
    }
    
    d.boos = d[tau,]
    
    e = d
    
    for (k in 1:(m-1)){
      
      e[,k] = d.boos[,k]-mean(d[,k])*(mean(d[,k]) >= - sqrt(w[k]^2/obs*2*log(log(obs))))
      
    }
    
    stat.boos[r] = max(0,max(sqrt(obs)*colMeans(e)/w))
  }
  
  p.value = mean((stat.boos > stat))
  
  list(p.value = p.value,stat.boos = stat.boos,stat = stat)
}

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# -----------------------------------------------
# BTC
# -----------------------------------------------

# Convert realized variance sequence into realized volatility sequence
btc_rv_df$realized_vol <- sqrt(btc_rv_df$RV_1min)

# ------------------------------
# Align data
# ------------------------------

# GARCH(1,1)
# Load volatility forecasts
btc_vol_forecast_garch <- read.csv("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub/data/rolling_forecasts_btc_daily.csv")

# Adjust data format and structure
btc_vol_forecast_garch <- btc_vol_forecast_garch[-nrow(btc_vol_forecast_garch), ]
btc_vol_forecast_garch <- btc_vol_forecast_garch$Sigma
btc_vol_forecast_garch <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_vol_forecast_garch
)

# Check length and specifications of data
head(btc_vol_forecast_garch)
tail(btc_vol_forecast_garch)

head(btc_rv_df)
tail(btc_rv_df)

length(btc_vol_forecast_garch$Date) == length(btc_rv_df$Date)
all(btc_vol_forecast_garch$Date == btc_rv_df$Date)


# EGARCH(,)
# Load volatility forecasts
btc_vol_forecast_egarch <- read.csv("")

# Adjust data format and structure

# Check length and specifications of data


# SGARCH(,)
# Load volatility forecasts
btc_vol_forecast_sgarch <- read.csv("")

# Adjust data format and structure

# Check length and specifications of data


# DCCGARCH(,)
# Load volatility forecasts
btc_vol_forecast_dccgarch <- read.csv("")

# Adjust data format and structure

# Check length and specifications of data


# Create list with all model forecasts
forecast_list_btc <- list(btc_vol_forecast_garch,
                          btc_vol_forecast_egarch,
                          btc_vol_forecast_sgarch,
                          btc_vol_forecast_dccgarch
                          )
# Define model names
model_names_btc <- c("GARCH(1,1)", "E GARCH(,)", "M GARCH(,)", "DCC GARCH (,)")

# ------------------------------
# Visual comparison
# ------------------------------

# Scatter plot comparison
for (i in 1:4){
  plot(forecast_list_btc[[i]]$forecast_vol, btc_rv_df$realized_vol,
      xlab = paste(model_names_btc[i], " predictions"), 
      ylab = "Volatility proxy")
  abline(0, 1, lty = 3)  
}

# Time-series plot overlap
for (i in 1:4){
  plot(btc_rv_df$Date, btc_rv_df$realized_vol, type = "l", col = "grey",
       xlab = "Date", ylab = "Volatility", main = paste("Volatility:", model_names_btc[i], " Forecast vs. Realized"))
  lines(forecast_list_btc[[i]]$Date, forecast_list_btc[[i]]$forecast_vol, col = "blue")
  
  legend("topright",
         legend = c("Realized Volatility", "Forecasted Volatility"),
         col = c("grey", "blue"), lty = 1, bty = "n",
         cex = 0.8, inset = c(0.02, 0.02))
}

# ------------------------------
# Mincer-Zarnowitz (MZ) Regression Test
# ------------------------------

# Classical MZ regression with robust White standard errors
regression_results_btc <- list()
robust_se_btc <- list()
coefs_btc <- list()
t_alpha_btc <- list()
t_beta_btc <- list()
p_alpha_btc <- list()
p_beta_btc <- list()

for (i in 1:4){
  
  # Run regression
  regression_results_btc[[i]] = lm(btc_rv_df$realized_vol ~ forecast_list_btc[[i]]$forecast_vol) 
  
  # Get robust (White) standard errors
  robust_se_btc[[i]] <- sqrt(diag(hccm(regression_results_btc[[i]], type = "hc0")))
  
  # Get coefficients
  coefs_btc[[i]] <- coef(regression_results_btc[[i]])
  
  # Compute t-stats for H0: alpha = 0 and H0: beta = 1
  t_alpha_btc[[i]] <- coefs_btc[[i]][1] / robust_se_btc[[i]][1]
  t_beta_btc[[i]]  <- (coefs_btc[[i]][2] - 1) / robust_se_btc[[i]][2]  # <-- test for beta = 1
  
  # Compute corresponding p-values
  p_alpha_btc[[i]] <- 2 * (1 - pt(abs(t_alpha_btc[[i]]), df = regression_results_btc[[i]]$df.residual))
  p_beta_btc[[i]]  <- 2 * (1 - pt(abs(t_beta_btc[[i]]), df = regression_results_btc[[i]]$df.residual))
  
  # Print
  cat(sprintf("α = %.4f (t = %.2f)%s, p = %.4f\n", coefs_btc[[i]][1], t_alpha_btc[[i]], ifelse(p_alpha_btc[[i]] < 0.01, " *", ""), p_alpha_btc[[i]]))
  cat(sprintf("β = %.4f (t for H0: β=1 = %.2f)%s, p = %.4f\n", coefs_btc[[i]][2], t_beta_btc[[i]], ifelse(p_beta_btc[[i]] < 0.01, " *", ""), p_beta_btc[[i]]))  
}

# ------------------------------
# Superior predictive ability (SPA)
# ------------------------------

# Data preparation
# Create u a vector with all different univariate GARCH predictions: each column, one model prediction
u <- cbind(garch11 = btc_vol_forecast_garch11$forecast_vol,
           e_garchXX = btc_vol_forecast_egarchXX$forecast_vol,
           m_garchXX = btc_vol_forecast_mgarchXX$forecast_vol,
           dcc_garchXX = btc_vol_forecast_dccgarchXX$forecast_vol
)

n <- dim(u)[1]

# Fill performance matrix with values of loss function
# Standard QLIKE loss function
perf_qlike <- matrix(,n,4)
for (j in 1:4)
{
  perf_qlike[, j] <- log(u[, j]) + btc_rv_df$realized_vol / u[, j] # Standard QLIKE
  #perf[,j] <- n*(-log(dnorm(us.stock.data$returns[, 6], mean = mean(us.stock.data$returns[, 6]), sd = sqrt(u[, j])))) # pseudo QLIKE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j]) #MAE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j])^2 #MSE
}

colMeans(perf_qlike)

# MSE loss function
perf_mse <- matrix(,n,4)
for (j in 1:4)
{
  #perf[, j] <- log(u[, j]) + btc_rv_df$realized_vol / u[, j] # Standard QLIKE
  #perf[,j] <- n*(-log(dnorm(us.stock.data$returns[, 6], mean = mean(us.stock.data$returns[, 6]), sd = sqrt(u[, j])))) # pseudo QLIKE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j]) #MAE
  perf_mse[, j] <- abs(btc_rv_df$realized_vol - u[, j])^2 #MSE
}

colMeans(perf_mse)

# Create list with both loss function matrices
perf <- list(perf_qlike, perf_mse)
loss_names <- c("QLIKE", "MSE")
model_names <- c("GARCH(1,1)", "E GARCH(,)", "M GARCH(,)", "DCC GARCH (,)")

# Run the SPA test, looping over each model as the benchmark using both QLIKE and MSE
for (loss_idx in 1:2){
  
  cat("\n===========================\n")
  cat("Loss function:", loss_names[loss_idx], "\n")
  cat("===========================\n")
  
  for (k in 1:4){
    
    cat("\nBenchmark model:", model_names[k], "\n")
    
    d = rep(0,1000)
    s = 1
    
    for (i in 1:10){ 
      
      #print(i)
      e = spa(per = perf[[loss_idx]], bench = k, m = 4, obs = n, q = 0.25, iter = 100, periodogram = T)
      d[s:(s+99)] = e$stat.boos
      s = s+100
    }
    
    p_value <- mean(d > e$stat)
    cat("SPA p-value:", round(p_value, 4), "\n")
    
  }
}

# ------------------------------
# Model confidence set (MCS)
# ------------------------------





# -----------------------------------------------
# ETH
# -----------------------------------------------

































