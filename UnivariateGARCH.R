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

ret <- na.omit(read.zoo("data/train_returns.csv", header = TRUE, sep = ",", format = "%Y-%m-%d", index.column = "Date"))
ret <- as.xts(ret)


# Simple example using an ARMA (1,0) GARCH (1,1) for the BTC return series
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)),
                    distribution.model = "std")
fit <- ugarchfit(spec = model, data = ret[, "BTC"])


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