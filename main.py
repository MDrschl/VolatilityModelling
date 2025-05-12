import pandas as pd
import numpy as np
import functionsCleaningDescriptives as dataFunc
import functionsUnivariate as funcUni

###########
# Data Loading, Preprocessing, Descriptive Statistics
###########
btc_full = dataFunc.load_close_series("data/BTCUSDT_1m.csv").sort_index()
eth_full = dataFunc.load_close_series("data/ETHUSDT_1m.csv").sort_index()

split_date = "2023-10-31"

# Train set: up to and including 2023-10-31
btc_train = btc_full.loc[:split_date]
eth_train = eth_full.loc[:split_date]

# Test set: after 2023-10-31
btc_test = btc_full.loc[split_date:]
eth_test = eth_full.loc[split_date:]

# Resample to daily prices and returns
btc_daily = btc_full["Close"].resample("1D").last().dropna()
eth_daily = eth_full["Close"].resample("1D").last().dropna()
btc_ret = btc_daily.pct_change().dropna()
eth_ret = eth_daily.pct_change().dropna()

# Daily prices
dataFunc.plot_series({"BTC": btc_daily, "ETH": eth_daily}, "Daily Close Price", "Price (USDT)")

# Daily returns
dataFunc.plot_series({"BTC": btc_ret, "ETH": eth_ret}, "Daily Returns", "Return")

freqs = {
    1440: "Daily",
    720: "12-hour",
    360: "6-hour",
    180: "3-hour",
    60: "1-hour"
}

results = []

for minutes, label in freqs.items():
    btc_resampled = btc_full["Close"].resample(f"{minutes}min").last().dropna()
    btc_ret_freq = btc_resampled.pct_change().dropna()
    btc_ret_freq.name = minutes
    btc_stats = dataFunc.compute_stats(btc_ret_freq)
    results.append(pd.Series(btc_stats, name=f"{label} return – BTC"))

for minutes, label in freqs.items():
    eth_resampled = eth_full["Close"].resample(f"{minutes}min").last().dropna()
    eth_ret_freq = eth_resampled.pct_change().dropna()
    eth_ret_freq.name = minutes
    eth_stats = dataFunc.compute_stats(eth_ret_freq)
    results.append(pd.Series(eth_stats, name=f"{label} return – ETH"))

summary_all = pd.DataFrame(results)
print(summary_all.round(3))

vol_btc = dataFunc.daily_ann_vol(btc_full["Close"].dropna())
vol_eth = dataFunc.daily_ann_vol(eth_full["Close"].dropna())

# Daily annualized volatility
dataFunc.plot_series({"BTC": vol_btc, "ETH": vol_eth}, "Daily Annualized Volatility (from 1-Minute Data)", "Volatility (%)")

# Autocorrelograms
dataFunc.plot_autocorrelogram(btc_ret, title="BTC Return Autocorrelogram")
dataFunc.plot_autocorrelogram(btc_ret, squared=True, title="BTC Squared Return Autocorrelogram")

dataFunc.plot_autocorrelogram(eth_ret, title="ETH Return Autocorrelogram")
dataFunc.plot_autocorrelogram(eth_ret, squared=True, title="ETH Squared Return Autocorrelogram")

# Diagnostics
all_returns = {}

for minutes, label in freqs.items():
    btc_resampled = btc_full["Close"].resample(f"{minutes}min").last().dropna()
    eth_resampled = eth_full["Close"].resample(f"{minutes}min").last().dropna()
    btc_ret_freq = btc_resampled.pct_change().dropna()
    eth_ret_freq = eth_resampled.pct_change().dropna()

    all_returns[f"{label} – BTC"] = btc_ret_freq
    all_returns[f"{label} – ETH"] = eth_ret_freq

diagnostics = dataFunc.test_distribution_diagnostics(all_returns)
print(diagnostics.round(3))

# Robust Diagnostics

# LM test

###########
# Univariate Models: In Sample
###########

# BTC
btc_returns_train = btc_train["Close"].resample("1D").last().pct_change().dropna()
funcUni.robust_acf_plot(btc_returns_train, title="BTC Return – Robust ACF")

btc_results, btc_best_model = funcUni.fit_garch_models(btc_returns_train, dist="t")
print("BTC GARCH Model Comparison:\n", btc_results.sort_values("AIC").round(3))

btc_params, btc_info = funcUni.summarize_garch_model(btc_best_model)
print("\nBTC Model Parameters:\n", btc_params)
print("\nBTC Fit Statistics:\n", btc_info)



# ETH
eth_returns_train = eth_train["Close"].resample("1D").last().pct_change().dropna()
funcUni.robust_acf_plot(eth_returns_train, title="ETH Return – Robust ACF")

eth_results, eth_best_model = funcUni.fit_garch_models(eth_returns_train, dist="t")
print("\nETH GARCH Model Comparison:\n", eth_results.sort_values("AIC").round(3))

eth_params, eth_info = funcUni.summarize_garch_model(eth_best_model)
print("\nETH Model Parameters:\n", eth_params)
print("\nETH Fit Statistics:\n", eth_info)

# Threshold GARCH
# Regime Switching

# Residual analysis: Robust ACF, Robust LB test

# Hourly, 6 hour frequency

###########
# Univariate Models: Out of Sample
###########

###########
# Multivariate Models: In Sample
###########

# VARMA-DCC-AGARCH
# Regime Switching

###########
# Multivariate Models: Out of Sample
###########