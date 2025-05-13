library(readr)
library(xts)
library(rmgarch)

# Load data
returns <- read_csv("train_returns.csv")
returns_xts <- xts(returns[, -1], order.by = as.Date(returns[[1]]))


# Define univariate GARCH specification (same for both series)
uspec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)

# Combine into a multivariate DCC-GARCH spec
spec <- dccspec(
  uspec = multispec(replicate(2, uspec)),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

# Fit the model
fit <- dccfit(spec, data = returns_xts)
