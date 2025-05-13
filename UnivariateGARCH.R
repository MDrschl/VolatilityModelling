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

btc_ret <- read.csv("data/train_returns.csv")
btc_ret$Date <- as.Date(btc_ret$Date)

ret <- na.omit(xts(btc_ret[, c("BTC", "ETH")], order.by = btc_ret$Date))

max_p <- 3
max_q <- 3
max_r <- 3
max_s <- 3

plot_list <- list()


# Initialize results tables
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
#par(mfrow = c(3,3))


for (col in colnames(ret)) {
  series <- ret[, col]
  
  #Initialise variables to store the best model's information
  best_model <- NULL
  best_aic <- Inf
  best_bic <- Inf
  best_ll <- -Inf
  best_order <- NULL
  aic_matrix <- matrix(NA, nrow = max_p, ncol = max_q)
  bic_matrix <- matrix(NA, nrow = max_p, ncol = max_q)
  llh_matrix <- matrix(NA, nrow = max_p, ncol = max_q)
  jb_val <- NA
  bp_val <- NA
  ml_val <- NA
  lb_val <- NA
  
  # Loop through different combinations of ARMA(p, q) and GARCH(r, s)
  for (p in 0:max_p) {
    for (q in 0:max_q) {
      for (r in 0:max_r) {
        for (s in 0:max_s) {
          if (p == 0 && q == 0 && r == 0 && s == 0) {
            next
          }
          
          # Fit ARMA-GARCH model
          model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(r, s)),
                              mean.model = list(armaOrder = c(p, q)),
                              distribution.model = "norm")
          fit <- tryCatch(ugarchfit(spec = model, data = series, solver = "nloptr"),
                          error = function(e) NULL)
          
          # Continue only if the model fitting was successful
          if (!is.null(fit)) {
            # Extract AIC, BIC, and calculate LL
            aic <- -2*fit@fit$LLH + 2*(r+s+p+q+1)
            bic <- -2*fit@fit$LLH + (r+s+p+q+1)*log(length(series))
            ll <- fit@fit$LLH
            
            aic_matrix[p, q] <- aic
            bic_matrix[p, q] <- bic
            llh_matrix[p, q] <- ll
            
            order <- paste("ARMA-GARCH(", p, ",", q, ";", r, ",", s, ")", sep = "")
            
            # Update the best model if necessary
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
  
  # Store the matrices for this series
  aic_matrices[[col]] <- aic_matrix
  bic_matrices[[col]] <- bic_matrix
  llh_matrices[[col]] <- llh_matrix
  
  # Store the best model for the current series
  best_models[[col]] <- best_model
  #stand_res <- best_model@fit$z
  #acf(stand_res^2, main = col, lag = 30) #For ACF of the Standardised squared residuals
  
  #plot_list[[col]] <- sigma_plot(series, best_model)
  
  # Append the results to the results_table
  results_table <- rbind(results_table, data.frame(
    Series = col,
    LLH = best_ll,
    AIC = best_aic,
    BIC = best_bic,
    Model = best_order
  ))
  # Append the tests to the test_table
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
coef(best_fit)
plot(best_fit, which = 2)
plot(best_fit, which = 3)

