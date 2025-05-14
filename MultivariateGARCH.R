library(rugarch)
library(zoo)
library(xts)

ret <- na.omit(read.zoo("data/train_returns.csv", header = TRUE, sep = ",", format = "%Y-%m-%d", index.column = "Date"))
ret <- as.xts(ret)
data <- ret[, c("BTC", "ETH")]

dcc_orders <- list(c(1, 0), c(0, 1), c(1, 1), c(2, 1), c(1, 2), c(2,2))
distributions <- c("mvnorm", "mvt")

results <- list()
best_aic <- Inf
best_model <- NULL
best_spec_desc <- ""

for (dcc_order in dcc_orders) {
  for (dist in distributions) {
    spec <- dccspec(
      uspec = multispec(list(
        ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(1, 0)),
                   distribution.model = "std"),
        ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(1, 0)),
                   distribution.model = "std"))),
      dccOrder = dcc_order,
      distribution = dist
    )
    fit <- tryCatch(dccfit(spec, data = data), error = function(e) NULL)
    
    if (!is.null(fit)) {
      aic_val <- infocriteria(fit)[1] # AIC is the first element
      spec_desc <- paste0("DCC(", dcc_order[1], ",", dcc_order[2], "), Dist: ", dist)
      results[[spec_desc]] <- list(fit = fit, aic = aic_val)
      
      if (aic_val < best_aic) {
        best_aic <- aic_val
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
print(lapply(results, function(res) if (is.list(res)) paste("AIC:", round(res$aic, 3)) else res))

print(paste("\nBest model based on AIC:", best_spec_desc, "with AIC:", round(best_aic, 3)))
coef(best_model)

correlation <- rcor(best_model)

# Extract dynamic correlation between BTC and ETH (off-diagonal)
btc_eth_corr <- correlation[1, 2, ]
btc_eth_corr <- xts(btc_eth_corr, order.by = index(ret))

autoplot(as.zoo(btc_eth_corr)) +
  labs(title = "Dynamic Correlation between BTC and ETH", y = "Correlation") +
  theme_minimal()

grangertest(ETH^2 ~ BTC^2, order = 1, data = as.data.frame(ret))

# Checking for vola spillover, but there might be more sophisticated checks for this
BTC_vol <- sigma(fit)[, 1]
ETH_vol <- sigma(fit)[, 2]

summary(lm(ETH_vol ~ lag(BTC_vol, 1)))

