library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom)

# Function to run test regressions for a given return series
run_bias_tests_on_returns <- function(ret, h = 1) {
  ret_lag <- lag(ret, h)
  
  df <- tibble(
    ret2 = ret^2,
    sign_bias = as.numeric(ret_lag < 0),
    neg_size_bias = ifelse(ret_lag < 0, ret_lag, 0),
    pos_size_bias = ifelse(ret_lag > 0, ret_lag, 0)
  ) %>% drop_na()
  
  list(
    sign = tidy(lm(ret2 ~ sign_bias, data = df))[2, ],
    neg  = tidy(lm(ret2 ~ neg_size_bias, data = df))[2, ],
    pos  = tidy(lm(ret2 ~ pos_size_bias, data = df))[2, ]
  )
}

# Define assets and frequencies
assets <- c("btc", "eth")
frequencies <- c("day", "6 hours")
lags <- c(1, 5)

# Run tests
results <- expand.grid(Asset = toupper(assets), Frequency = frequencies, Lag = lags, stringsAsFactors = FALSE) %>%
  mutate(across(everything(), as.character))

# Add results columns
results <- results %>%
  mutate(
    Sign_Bias = NA_character_,
    Neg_Size_Bias = NA_character_,
    Pos_Size_Bias = NA_character_
  )

# Fill in table
for (i in 1:nrow(results)) {
  asset <- tolower(results$Asset[i])
  freq <- results$Frequency[i]
  lag <- as.numeric(results$Lag[i])
  
  ret_df <- returns_model_fitting[[paste0(asset, "_train")]][[freq]]
  ret_xts <- xts(ret_df$Return, order.by = ret_df$Time)
  ret_num <- as.numeric(ret_xts)
  
  out <- run_bias_tests_on_returns(ret_num, h = lag)
  
  results$Sign_Bias[i] <- sprintf("%.3f%s", out$sign$estimate,
                                  ifelse(out$sign$p.value < 0.05, "*", ""))
  results$Neg_Size_Bias[i] <- sprintf("%.3f%s", out$neg$estimate,
                                      ifelse(out$neg$p.value < 0.05, "*", ""))
  results$Pos_Size_Bias[i] <- sprintf("%.3f%s", out$pos$estimate,
                                      ifelse(out$pos$p.value < 0.05, "*", ""))
}

# Pivot for nicer output like in slide
results_final <- results %>%
  pivot_longer(cols = c(Sign_Bias, Neg_Size_Bias, Pos_Size_Bias), names_to = "Test", values_to = "Value") %>%
  pivot_wider(names_from = c(Test, Lag), values_from = Value) %>%
  arrange(Asset, Frequency)


print(results_wide, row.names = FALSE)
