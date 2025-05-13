library(readr)
library(xts)
library(rmgarch)

# Load data
returns <- read_csv("train_returns.csv")
returns_xts <- xts(returns[, -1], order.by = as.Date(returns[[1]]))
