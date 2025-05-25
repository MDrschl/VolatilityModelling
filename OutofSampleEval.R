#-----------
# Set working directory
#-----------
setwd("/Users/nathanielsuchin/Library/Mobile Documents/com~apple~CloudDocs/Documents/University/University St. Gallen/2025 Spring Semester/Financial Volatility/Group Assignment/GitHub") # replace with own path

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
library(MCS)

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
      
      # WITH THIS (corrected):
      if (mean(d[,k]) >= -sqrt(w[k]^2/obs*2*log(log(obs)))) {
        # Standard recentering when condition is met
        e[,k] = d.boos[,k] - mean(d[,k])
      } else {
        # Modified recentering when benchmark is very superior
        e[,k] = d.boos[,k] - mean(d.boos[,k])
      }
      
    }
    
    stat.boos[r] = max(0,max(sqrt(obs)*colMeans(e)/w))
  }
  
  p.value = mean((stat.boos > stat))
  
  list(p.value = p.value,stat.boos = stat.boos,stat = stat)
}

#-----------
# MCF test
#-----------
# CORRECTED MCS FUNCTION - REPLACE YOUR EXISTING ONE!
mcs_test <- function(losses, alpha = 0.05, B = 1000, model_names = NULL) {
  # losses: matrix where each column is a model's loss series
  # alpha: significance level  
  # B: number of bootstrap replications
  # model_names: optional vector of model names for display
  
  n <- nrow(losses)
  k <- ncol(losses)
  
  # Use default names if not provided
  if (is.null(model_names)) {
    model_names <- paste("Model", 1:k)
  }
  
  # Initialize
  M <- 1:k  # Set of models still in contention
  eliminated <- c()
  
  cat("Starting MCS with", k, "models\n")
  
  while(length(M) > 1) {
    
    cat("Testing", length(M), "models:", M, "\n")
    
    # Calculate pairwise loss differences for models in M
    if(length(M) == 2) {
      d_matrix <- matrix(losses[, M[1]] - losses[, M[2]], ncol = 1)
      model_pairs <- matrix(c(M[1], M[2]), nrow = 1)
    } else {
      d_matrix <- matrix(0, n, choose(length(M), 2))
      col_idx <- 1
      model_pairs <- matrix(0, choose(length(M), 2), 2)
      
      for(i in 1:(length(M)-1)) {
        for(j in (i+1):length(M)) {
          d_matrix[, col_idx] <- losses[, M[i]] - losses[, M[j]]
          model_pairs[col_idx, ] <- c(M[i], M[j])
          col_idx <- col_idx + 1
        }
      }
    }
    
    # Calculate test statistics
    d_bar <- colMeans(d_matrix)
    
    # Bootstrap variance estimation
    bootstrap_vars <- numeric(ncol(d_matrix))
    for(b in 1:B) {
      boot_indices <- sample(1:n, n, replace = TRUE)
      boot_d <- d_matrix[boot_indices, , drop = FALSE]
      boot_d_bar <- colMeans(boot_d)
      bootstrap_vars <- bootstrap_vars + (boot_d_bar - d_bar)^2
    }
    bootstrap_vars <- bootstrap_vars / B
    
    # Calculate standardized test statistics
    t_stats <- abs(d_bar) / sqrt(bootstrap_vars + 1e-8)
    
    # Max standardized statistic
    t_max <- max(t_stats)
    
    # Bootstrap p-values
    bootstrap_t_max <- numeric(B)
    for(b in 1:B) {
      boot_indices <- sample(1:n, n, replace = TRUE)
      boot_d <- d_matrix[boot_indices, , drop = FALSE]
      boot_d_bar <- colMeans(boot_d)
      boot_t_stats <- abs(boot_d_bar - d_bar) / sqrt(bootstrap_vars + 1e-8)
      bootstrap_t_max[b] <- max(boot_t_stats)
    }
    
    # P-value
    p_value <- mean(bootstrap_t_max > t_max)
    
    cat("Max test statistic:", round(t_max, 4), "P-value:", round(p_value, 4), "\n")
    
    if(p_value > alpha) {
      # Cannot reject that all remaining models are equivalent
      cat("Cannot reject equivalence of remaining models\n")
      break
    } else {
      # Eliminate worst performing model
      worst_idx <- which.max(abs(d_bar))
      worst_pair <- model_pairs[worst_idx, ]
      
      # Determine which model in the pair performs worse (higher loss is worse)
      model_i_avg_loss <- mean(losses[, worst_pair[1]])
      model_j_avg_loss <- mean(losses[, worst_pair[2]])
      
      if(model_i_avg_loss > model_j_avg_loss) {
        eliminate <- worst_pair[1]
      } else {
        eliminate <- worst_pair[2]
      }
      
      # FIXED: Use the provided model_names parameter instead of hardcoded variable
      cat("Eliminating model", eliminate, "(", model_names[eliminate], ")\n")
      eliminated <- c(eliminated, eliminate)
      M <- setdiff(M, eliminate)
    }
  }
  
  list(
    MCS = M,
    eliminated = eliminated,
    eliminated_order = eliminated
  )
}

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

# -----------------------------------------------
# BTC
# -----------------------------------------------

# Load BTC realized variances
btc_rv_df <- read.csv("btc_rv_df.csv")
btc_rv_df$Date <- as.Date(btc_rv_df$Date)

# Convert realized variance sequence into realized volatility sequence
btc_rv_df$realized_vol <- sqrt(btc_rv_df$RV_1min)

# ------------------------------
# Align data
# ------------------------------

# GARCH(1,1)
# Load volatility forecasts
btc_vol_forecast_garch <- read.csv("data/Forecasts/btc_garch_forecasts.csv")

# Adjust data format and structure
btc_vol_forecast_garch <- btc_vol_forecast_garch[-nrow(btc_vol_forecast_garch), ]
btc_vol_forecast_garch <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_vol_forecast_garch[[2]]
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
btc_vol_forecast_egarch <- read.csv("data/Forecasts/btc_egarch_forecasts.csv")

# Adjust data format and structure
btc_vol_forecast_egarch <- btc_vol_forecast_egarch[-nrow(btc_vol_forecast_egarch), ]
btc_vol_forecast_egarch <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_vol_forecast_egarch[[2]]
)

# Check length and specifications of data
length(btc_vol_forecast_egarch$Date) == length(btc_rv_df$Date)
all(btc_vol_forecast_egarch$Date == btc_rv_df$Date)

# MSGARCH(,)
# Load volatility forecasts
btc_vol_forecast_msgarch <- read.csv("data/Forecasts/btc_msgarch_forecasts.csv")

# Adjust data format and structure
btc_vol_forecast_msgarch <- btc_vol_forecast_msgarch[-nrow(btc_vol_forecast_msgarch), ]
btc_vol_forecast_msgarch <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_vol_forecast_msgarch[[2]]
)

# Check length and specifications of data
length(btc_vol_forecast_msgarch$Date) == length(btc_rv_df$Date)
all(btc_vol_forecast_msgarch$Date == btc_rv_df$Date)

# DCC-GARCH(,)
# Load volatility forecasts
btc_vol_forecast_dccgarch <- read.csv("data/Forecasts/btc_eth_dcc_forecast.csv")

# Adjust data format and structure
btc_vol_forecast_dccgarch <- btc_vol_forecast_dccgarch[-nrow(btc_vol_forecast_dccgarch), ]
btc_vol_forecast_dccgarch <- data.frame(
  Date = btc_rv_df$Date,
  forecast_vol = btc_vol_forecast_dccgarch$BTC_StdDev
)

# Check length and specifications of data
length(btc_vol_forecast_dccgarch$Date) == length(btc_rv_df$Date)
all(btc_vol_forecast_dccgarch$Date == btc_rv_df$Date)

# Create list with all model forecasts
forecast_list_btc <- list(btc_vol_forecast_garch,
                          btc_vol_forecast_egarch,
                          btc_vol_forecast_msgarch,
                          btc_vol_forecast_dccgarch
                          )
#for (i in 1:4) {
#  forecast_list_btc[[i]]$Date <- as.Date(forecast_list_btc[[i]]$Date)
#}

# Define model names
model_names_btc <- c("GARCH(1,1)", "EGARCH(1,1)", "MSGARCH(1,1)", "DCC-GARCH(1,1)")

# ------------------------------
# Visual comparison
# ------------------------------

# Scatter plot comparison

for (i in 1:4){
  plot(forecast_list_btc[[i]]$forecast_vol, btc_rv_df$realized_vol,
      xlab = paste(model_names_btc[i], " predictions"), 
      ylab = "Volatility proxy",
      main = "BTC")
  abline(0, 1, lty = 3)  
}



# Set up the plotting area for 1 row and 4 columns
par(mfrow=c(1, 4))

# Your plotting loop
for (i in 1:4){
  plot(forecast_list_btc[[i]]$forecast_vol, btc_rv_df$realized_vol,
       xlab = paste(model_names_btc[i], " predictions"), 
       ylab = "Volatility proxy",
       main = "BTC")
  abline(0, 1, lty = 3)  
}

# Reset to default single plot
par(mfrow=c(1, 1))

# Time-series plot overlap
for (i in 1:4){
  plot(btc_rv_df$Date, btc_rv_df$realized_vol, type = "l", col = "grey",
       xlab = "Date", ylab = "Volatility", main = paste("Volatility:", model_names_btc[i], " Forecast vs. Realized (BTC)"))
  lines(forecast_list_btc[[i]]$Date, forecast_list_btc[[i]]$forecast_vol, col = "blue")
  
  legend("topright",
         legend = c("Realized Volatility", "Forecasted Volatility"),
         col = c("grey", "blue"), lty = 1, bty = "n",
         cex = 0.6, inset = c(0.01, 0.01))
}



# Solution 1: Reduce margins significantly
par(mfrow=c(4, 1), mar=c(2, 4, 2, 1))  # smaller margins: bottom, left, top, right

# Time-series plot overlap
for (i in 1:4){
  plot(btc_rv_df$Date, btc_rv_df$realized_vol, type = "l", col = "grey",
       xlab = "Date", ylab = "Volatility", main = paste("Volatility:", model_names_btc[i], " Forecast vs. Realized (BTC)"))
  lines(forecast_list_btc[[i]]$Date, forecast_list_btc[[i]]$forecast_vol, col = "blue")
  
  legend("topright",
         legend = c("Realized Volatility", "Forecasted Volatility"),
         col = c("grey", "blue"), lty = 1, bty = "n",
         cex = 0.6, inset = c(0.01, 0.01))
}

par(mfrow=c(1, 1))




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

}

# Summary table
# Initialize matrix to store results
summary_table_mz <- matrix(NA, nrow = 4, ncol = 4)
rownames(summary_table_mz) <- c("alpha (p-value)", "alpha t-stat", "beta (p-value)", "beta t-stat")
colnames(summary_table_mz) <- c("GARCH(1,1)", "EGARCH(1,1)", "MSGARCH(1,1)", "DCC-GARCH(1,1)")

# Fill in the matrix
for (i in 1:4) {
  alpha_val <- round(coefs_btc[[i]][1], 4)
  alpha_p   <- round(p_alpha_btc[[i]], 4)
  beta_val  <- round(coefs_btc[[i]][2], 4)
  beta_p    <- round(p_beta_btc[[i]], 4)
  
  summary_table_mz[1, i] <- paste0(alpha_val, " (", alpha_p, ")")
  summary_table_mz[2, i] <- round(t_alpha_btc[[i]], 2)
  summary_table_mz[3, i] <- paste0(beta_val, " (", beta_p, ")")
  summary_table_mz[4, i] <- round(t_beta_btc[[i]], 2)
}

# Convert to data frame for cleaner display
summary_table_mz <- as.data.frame(summary_table_mz)

# Print
print(summary_table_mz)

# ------------------------------
# Superior predictive ability (SPA)
# ------------------------------

# Data preparation
# Create u a vector with all different univariate GARCH predictions: each column, one model prediction
u <- cbind(
  garch11     = forecast_list_btc[[1]]$forecast_vol,
  egarch      = forecast_list_btc[[2]]$forecast_vol,
  msgarch     = forecast_list_btc[[3]]$forecast_vol
  #dcc_garch   = forecast_list_btc[[4]]$forecast_vol
)

n <- dim(u)[1]

# Fill performance matrix with values of loss function
# Standard QLIKE loss function
perf_qlike <- matrix(,n,3)
for (j in 1:3)
{
  perf_qlike[, j] <- -log(u[, j]) + btc_rv_df$realized_vol / u[, j] # Standard QLIKE
  #perf[,j] <- n*(-log(dnorm(us.stock.data$returns[, 6], mean = mean(us.stock.data$returns[, 6]), sd = sqrt(u[, j])))) # pseudo QLIKE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j]) #MAE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j])^2 #MSE
}

colMeans(perf_qlike)

# MSE loss function
perf_mse <- matrix(,n,3)
for (j in 1:3)
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

# Run the SPA test, looping over each model as the benchmark using both QLIKE and MSE
for (loss_idx in 1:2){
  
  cat("\n===========================\n")
  cat("Loss function:", loss_names[loss_idx], "\n")
  cat("===========================\n")
  
  for (k in 1:3){
    cat("\nBenchmark model:", model_names_btc[k], "\n")
    
    d = rep(0, 10000)  # 1000 * 10 = 10000 total
    #          ^^^^^ Increase this
    s = 1
    
    for (i in 1:10){ 
      e = spa(per = perf[[loss_idx]], bench = k, m = 3, obs = n, q = 0.25, iter = 1000, periodogram = T)
      #                                                                      ^^^^ Keep at 1000
      d[s:(s+999)] = e$stat.boos   # Extract all 1000 values
      #     ^^^ Change to 999
      s = s + 1000  # Increment by 1000
      #       ^^^^ Change to 1000
    }
    
    p_value <- mean(d > e$stat)
    cat("SPA p-value:", round(p_value, 4), "\n")
  }
}

# Calculate mean loss values
qlike_means_btc <- colMeans(perf_qlike)
mse_means_btc <- colMeans(perf_mse)

# Create a summary table with full precision for MSE
loss_summary_btc <- data.frame(
  Model = model_names_btc[1:3],  # Only first 3 models since DCC is commented out
  QLIKE = round(qlike_means_btc, 4),
  MSE = round(mse_means_btc, 6)  # 6 decimal places for MSE
)

# Display the table
print(loss_summary_btc)

# ------------------------------
# Model confidence set (MCS)
# ------------------------------

# MCS for QLIKE
cat("--- MCS for QLIKE Loss ---\n")
mcs_qlike_manual <- mcs_test(perf_qlike, alpha = 0.05, B = 500, model_names = model_names_btc[1:3])
cat("\nMCS for QLIKE (5% level):\n")
cat("Models in MCS:", paste(model_names_btc[1:3][mcs_qlike_manual$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_btc[1:3][mcs_qlike_manual$eliminated], collapse = ", "), "\n")

cat("\n")
mcs_qlike_manual_10 <- mcs_test(perf_qlike, alpha = 0.10, B = 500, model_names = model_names_btc[1:3])
cat("\nMCS for QLIKE (10% level):\n")
cat("Models in MCS:", paste(model_names_btc[1:3][mcs_qlike_manual_10$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_btc[1:3][mcs_qlike_manual_10$eliminated], collapse = ", "), "\n")

# MCS for MSE
cat("\n--- MCS for MSE Loss ---\n")
mcs_mse_manual <- mcs_test(perf_mse, alpha = 0.05, B = 500, model_names = model_names_btc[1:3])
cat("\nMCS for MSE (5% level):\n")
cat("Models in MCS:", paste(model_names_btc[1:3][mcs_mse_manual$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_btc[1:3][mcs_mse_manual$eliminated], collapse = ", "), "\n")

cat("\n")
mcs_mse_manual_10 <- mcs_test(perf_mse, alpha = 0.10, B = 500, model_names = model_names_btc[1:3])
cat("\nMCS for MSE (10% level):\n")
cat("Models in MCS:", paste(model_names_btc[1:3][mcs_mse_manual_10$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_btc[1:3][mcs_mse_manual_10$eliminated], collapse = ", "), "\n")

# -----------------------------------------------
# ETH
# -----------------------------------------------

# Load ETH realized variances
eth_rv_df <- read.csv("eth_rv_df.csv")
eth_rv_df$Date <- as.Date(eth_rv_df$Date)

# Convert realized variance sequence into realized volatility sequence
eth_rv_df$realized_vol <- sqrt(eth_rv_df$RV_1min)

# ------------------------------
# Align data
# ------------------------------

# GARCH(1,1)
# Load volatility forecasts
eth_vol_forecast_garch <- read.csv("data/Forecasts/eth_garch_forecasts.csv")

# Adjust data format and structure
eth_vol_forecast_garch <- eth_vol_forecast_garch[-nrow(eth_vol_forecast_garch), ]
eth_vol_forecast_garch <- data.frame(
  Date = eth_rv_df$Date,
  forecast_vol = eth_vol_forecast_garch[[2]]
)

# Check length and specifications of data
head(eth_vol_forecast_garch)
tail(eth_vol_forecast_garch)

head(eth_rv_df)
tail(eth_rv_df)

length(eth_vol_forecast_garch$Date) == length(eth_rv_df$Date)
all(eth_vol_forecast_garch$Date == eth_rv_df$Date)

# EGARCH(1,1)
# Load volatility forecasts
eth_vol_forecast_egarch <- read.csv("data/Forecasts/eth_egarch_forecasts.csv")

# Adjust data format and structure
eth_vol_forecast_egarch <- eth_vol_forecast_egarch[-nrow(eth_vol_forecast_egarch), ]
eth_vol_forecast_egarch <- data.frame(
  Date = eth_rv_df$Date,
  forecast_vol = eth_vol_forecast_egarch[[2]]
)

# Check length and specifications of data
length(eth_vol_forecast_egarch$Date) == length(eth_rv_df$Date)
all(eth_vol_forecast_egarch$Date == eth_rv_df$Date)

# MSGARCH(1,1)
# Load volatility forecasts
eth_vol_forecast_msgarch <- read.csv("data/Forecasts/eth_msgarch_forecasts.csv")

# Adjust data format and structure
eth_vol_forecast_msgarch <- eth_vol_forecast_msgarch[-nrow(eth_vol_forecast_msgarch), ]
eth_vol_forecast_msgarch <- data.frame(
  Date = eth_rv_df$Date,
  forecast_vol = eth_vol_forecast_msgarch[[2]]
)

# Check length and specifications of data
length(eth_vol_forecast_msgarch$Date) == length(eth_rv_df$Date)
all(eth_vol_forecast_msgarch$Date == eth_rv_df$Date)

# DCC-GARCH(1,1)
# Load volatility forecasts
eth_vol_forecast_dccgarch <- read.csv("data/Forecasts/btc_eth_dcc_forecast.csv")

# Adjust data format and structure
eth_vol_forecast_dccgarch <- eth_vol_forecast_dccgarch[-nrow(eth_vol_forecast_dccgarch), ]
eth_vol_forecast_dccgarch <- data.frame(
  Date = eth_rv_df$Date,
  forecast_vol = eth_vol_forecast_dccgarch$ETH_StdDev
)

# Check length and specifications of data
length(eth_vol_forecast_dccgarch$Date) == length(eth_rv_df$Date)
all(eth_vol_forecast_dccgarch$Date == eth_rv_df$Date)

# Create list with all model forecasts
forecast_list_eth <- list(eth_vol_forecast_garch,
                          eth_vol_forecast_egarch,
                          eth_vol_forecast_msgarch,
                          eth_vol_forecast_dccgarch
)

# Define model names
model_names_eth <- c("GARCH(1,1)", "EGARCH(1,1)", "MSGARCH(1,1)", "DCC-GARCH(1,1)")

# ------------------------------
# Visual comparison
# ------------------------------

# Scatter plot comparison
for (i in 1:4){
  plot(forecast_list_eth[[i]]$forecast_vol, eth_rv_df$realized_vol,
       xlab = paste(model_names_eth[i], " predictions"), 
       ylab = "Volatility proxy",
       main = "ETH")
  abline(0, 1, lty = 3)  
}

# Scatter plot comparison - 4 plots side by side
par(mfrow=c(1, 4), mar=c(5, 4, 2, 1), cex=0.7)  # Increased bottom margin, smaller text

for (i in 1:4){
  plot(forecast_list_eth[[i]]$forecast_vol, eth_rv_df$realized_vol,
       xlab = paste(model_names_eth[i], "predictions"), 
       ylab = "Volatility proxy",
       main = "ETH")
  abline(0, 1, lty = 3)  
}
par(mfrow=c(1, 1))  # Reset to default

# Time-series plot overlap
for (i in 1:4){
  plot(eth_rv_df$Date, eth_rv_df$realized_vol, type = "l", col = "grey",
       xlab = "Date", ylab = "Volatility", main = paste("Volatility:", model_names_eth[i], " Forecast vs. Realized (ETH)"))
  lines(forecast_list_eth[[i]]$Date, forecast_list_eth[[i]]$forecast_vol, col = "blue")
  
  legend("topright",
         legend = c("Realized Volatility", "Forecasted Volatility"),
         col = c("grey", "blue"), lty = 1, bty = "n",
         cex = 0.6, inset = c(0.01, 0.01))
}

# Time-series plot overlap - 4 plots stacked vertically
par(mfrow=c(4, 1), mar=c(2, 4, 2, 1))  # Reduced margins to avoid "margins too large" error

for (i in 1:4){
  plot(eth_rv_df$Date, eth_rv_df$realized_vol, type = "l", col = "grey",
       xlab = "Date", ylab = "Volatility", main = paste("Volatility:", model_names_eth[i], " Forecast vs. Realized (ETH)"))
  lines(forecast_list_eth[[i]]$Date, forecast_list_eth[[i]]$forecast_vol, col = "blue")
  
  legend("topright",
         legend = c("Realized Volatility", "Forecasted Volatility"),
         col = c("grey", "blue"), lty = 1, bty = "n",
         cex = 0.6, inset = c(0.01, 0.01))
}

par(mfrow=c(1, 1))  # Reset to default


# ------------------------------
# Mincer-Zarnowitz (MZ) Regression Test
# ------------------------------

# Classical MZ regression with robust White standard errors
regression_results_eth <- list()
robust_se_eth <- list()
coefs_eth <- list()
t_alpha_eth <- list()
t_beta_eth <- list()
p_alpha_eth <- list()
p_beta_eth <- list()

for (i in 1:4){
  
  # Run regression
  regression_results_eth[[i]] = lm(eth_rv_df$realized_vol ~ forecast_list_eth[[i]]$forecast_vol) 
  
  # Get robust (White) standard errors
  robust_se_eth[[i]] <- sqrt(diag(hccm(regression_results_eth[[i]], type = "hc0")))
  
  # Get coefficients
  coefs_eth[[i]] <- coef(regression_results_eth[[i]])
  
  # Compute t-stats for H0: alpha = 0 and H0: beta = 1
  t_alpha_eth[[i]] <- coefs_eth[[i]][1] / robust_se_eth[[i]][1]
  t_beta_eth[[i]]  <- (coefs_eth[[i]][2] - 1) / robust_se_eth[[i]][2]  # <-- test for beta = 1
  
  # Compute corresponding p-values
  p_alpha_eth[[i]] <- 2 * (1 - pt(abs(t_alpha_eth[[i]]), df = regression_results_eth[[i]]$df.residual))
  p_beta_eth[[i]]  <- 2 * (1 - pt(abs(t_beta_eth[[i]]), df = regression_results_eth[[i]]$df.residual))
  
}

# Summary table
# Initialize matrix to store results
summary_table_mz_eth <- matrix(NA, nrow = 4, ncol = 4)
rownames(summary_table_mz_eth) <- c("alpha (p-value)", "alpha t-stat", "beta (p-value)", "beta t-stat")
colnames(summary_table_mz_eth) <- c("GARCH(1,1)", "EGARCH(1,1)", "MSGARCH(1,1)", "DCC-GARCH(1,1)")

# Fill in the matrix
for (i in 1:4) {
  alpha_val <- round(coefs_eth[[i]][1], 4)
  alpha_p   <- round(p_alpha_eth[[i]], 4)
  beta_val  <- round(coefs_eth[[i]][2], 4)
  beta_p    <- round(p_beta_eth[[i]], 4)
  
  summary_table_mz_eth[1, i] <- paste0(alpha_val, " (", alpha_p, ")")
  summary_table_mz_eth[2, i] <- round(t_alpha_eth[[i]], 2)
  summary_table_mz_eth[3, i] <- paste0(beta_val, " (", beta_p, ")")
  summary_table_mz_eth[4, i] <- round(t_beta_eth[[i]], 2)
}

# Convert to data frame for cleaner display
summary_table_mz_eth <- as.data.frame(summary_table_mz_eth)

# Print
print(summary_table_mz_eth)

# ------------------------------
# Superior predictive ability (SPA)
# ------------------------------

# Data preparation
# Create u a vector with all different univariate GARCH predictions: each column, one model prediction
u_eth <- cbind(
  garch11     = forecast_list_eth[[1]]$forecast_vol,
  egarch      = forecast_list_eth[[2]]$forecast_vol,
  msgarch     = forecast_list_eth[[3]]$forecast_vol
  #dcc_garch   = forecast_list_eth[[4]]$forecast_vol
)

n_eth <- dim(u_eth)[1]

# Fill performance matrix with values of loss function
# Standard QLIKE loss function
perf_qlike_eth <- matrix(,n_eth,3)
for (j in 1:3)
{
  perf_qlike_eth[, j] <- -log(u_eth[, j]) + eth_rv_df$realized_vol / u_eth[, j] # Standard QLIKE
  #perf[,j] <- n*(-log(dnorm(us.stock.data$returns[, 6], mean = mean(us.stock.data$returns[, 6]), sd = sqrt(u[, j])))) # pseudo QLIKE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j]) #MAE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j])^2 #MSE
}

colMeans(perf_qlike_eth)

# MSE loss function
perf_mse_eth <- matrix(,n_eth,3)
for (j in 1:3)
{
  #perf[, j] <- log(u[, j]) + eth_rv_df$realized_vol / u[, j] # Standard QLIKE
  #perf[,j] <- n*(-log(dnorm(us.stock.data$returns[, 6], mean = mean(us.stock.data$returns[, 6]), sd = sqrt(u[, j])))) # pseudo QLIKE
  #perf[, j] <- abs(us.stock.data$real.cov[, 6, 6]-u[, j]) #MAE
  perf_mse_eth[, j] <- abs(eth_rv_df$realized_vol - u_eth[, j])^2 #MSE
}

colMeans(perf_mse_eth)

# Create list with both loss function matrices
perf_eth <- list(perf_qlike_eth, perf_mse_eth)
loss_names <- c("QLIKE", "MSE")

# Run the SPA test, looping over each model as the benchmark using both QLIKE and MSE
for (loss_idx in 1:2){
  
  cat("\n===========================\n")
  cat("Loss function:", loss_names[loss_idx], "\n")
  cat("===========================\n")
  
  for (k in 1:3){
    cat("\nBenchmark model:", model_names_eth[k], "\n")
    
    d = rep(0, 10000)  # 1000 * 10 = 10000 total
    #          ^^^^^ Increase this
    s = 1
    
    for (i in 1:10){ 
      e = spa(per = perf_eth[[loss_idx]], bench = k, m = 3, obs = n_eth, q = 0.25, iter = 1000, periodogram = T)
      #                                                                      ^^^^ Keep at 1000
      d[s:(s+999)] = e$stat.boos   # Extract all 1000 values
      #     ^^^ Change to 999
      s = s + 1000  # Increment by 1000
      #       ^^^^ Change to 1000
    }
    
    p_value <- mean(d > e$stat)
    cat("SPA p-value:", round(p_value, 4), "\n")
  }
}

# Calculate mean loss values
qlike_means_eth <- colMeans(perf_qlike_eth)
mse_means_eth <- colMeans(perf_mse_eth)

# Create a summary table with more decimal places for MSE
loss_summary_eth <- data.frame(
  Model = model_names_eth[1:3],
  QLIKE = round(qlike_means_eth, 4),
  MSE = mse_means_eth  # Don't round MSE
)

# Display with full precision
print(loss_summary_eth)

# Or if you want specific decimal places for MSE
loss_summary_eth <- data.frame(
  Model = model_names_eth[1:3],
  QLIKE = round(qlike_means_eth, 4),
  MSE = round(mse_means_eth, 6)  # 6 decimal places
)

print(loss_summary_eth)

# ------------------------------
# Model confidence set (MCS)
# ------------------------------

# For ETH:
# MCS for QLIKE
cat("--- MCS for QLIKE Loss ---\n")
mcs_qlike_manual_eth <- mcs_test(perf_qlike_eth, alpha = 0.05, B = 500, model_names = model_names_eth[1:3])
cat("\nMCS for QLIKE (5% level):\n")
cat("Models in MCS:", paste(model_names_eth[1:3][mcs_qlike_manual_eth$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_eth[1:3][mcs_qlike_manual_eth$eliminated], collapse = ", "), "\n")

cat("\n")
mcs_qlike_manual_10_eth <- mcs_test(perf_qlike_eth, alpha = 0.10, B = 500, model_names = model_names_eth[1:3])
cat("\nMCS for QLIKE (10% level):\n")
cat("Models in MCS:", paste(model_names_eth[1:3][mcs_qlike_manual_10_eth$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_eth[1:3][mcs_qlike_manual_10_eth$eliminated], collapse = ", "), "\n")

# MCS for MSE
cat("\n--- MCS for MSE Loss ---\n")
mcs_mse_manual_eth <- mcs_test(perf_mse_eth, alpha = 0.05, B = 500, model_names = model_names_eth[1:3])
cat("\nMCS for MSE (5% level):\n")
cat("Models in MCS:", paste(model_names_eth[1:3][mcs_mse_manual_eth$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_eth[1:3][mcs_mse_manual_eth$eliminated], collapse = ", "), "\n")

cat("\n")
mcs_mse_manual_10_eth <- mcs_test(perf_mse_eth, alpha = 0.10, B = 500, model_names = model_names_eth[1:3])
cat("\nMCS for MSE (10% level):\n")
cat("Models in MCS:", paste(model_names_eth[1:3][mcs_mse_manual_10_eth$MCS], collapse = ", "), "\n")
cat("Models eliminated:", paste(model_names_eth[1:3][mcs_mse_manual_10_eth$eliminated], collapse = ", "), "\n")
