import pandas as pd
import numpy as np
import functionsCleaningDescriptives as dataFunc
import functionsResidualsAnalysis as funcRes
import functionsUnivariate as funcUni

###########
# Data Loading, Preprocessing, Descriptive Statistics
###########

btc = dataFunc.load_close_series("data/BTCUSDT_1m.csv")
eth = dataFunc.load_close_series("data/ETHUSDT_1m.csv")

# Resample to daily prices and returns
btc_daily = btc["Close"].resample("1D").last().dropna()
eth_daily = eth["Close"].resample("1D").last().dropna()
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
    btc_resampled = btc["Close"].resample(f"{minutes}min").last().dropna()
    btc_ret_freq = btc_resampled.pct_change().dropna()
    btc_ret_freq.name = minutes
    btc_stats = dataFunc.compute_stats(btc_ret_freq)
    results.append(pd.Series(btc_stats, name=f"{label} return – BTC"))

for minutes, label in freqs.items():
    eth_resampled = eth["Close"].resample(f"{minutes}min").last().dropna()
    eth_ret_freq = eth_resampled.pct_change().dropna()
    eth_ret_freq.name = minutes
    eth_stats = dataFunc.compute_stats(eth_ret_freq)
    results.append(pd.Series(eth_stats, name=f"{label} return – ETH"))

summary_all = pd.DataFrame(results)
print(summary_all.round(3))

vol_btc = dataFunc.daily_ann_vol(btc["Close"].dropna())
vol_eth = dataFunc.daily_ann_vol(eth["Close"].dropna())

# Daily annualized volatility
dataFunc.plot_series({"BTC": vol_btc, "ETH": vol_eth}, "Daily Annualized Volatility (from 1-Minute Data)", "Volatility (%)")

# Autocorrelograms
dataFunc.plot_autocorrelogram(btc_ret, title="BTC Return Autocorrelogram")
dataFunc.plot_autocorrelogram(btc_ret, squared=True, title="BTC Squared Return Autocorrelogram")

dataFunc.plot_autocorrelogram(eth_ret, title="ETH Return Autocorrelogram")
dataFunc.plot_autocorrelogram(eth_ret, squared=True, title="ETH Squared Return Autocorrelogram")

all_returns = {}

for minutes, label in freqs.items():
    btc_resampled = btc["Close"].resample(f"{minutes}min").last().dropna()
    eth_resampled = eth["Close"].resample(f"{minutes}min").last().dropna()
    btc_ret_freq = btc_resampled.pct_change().dropna()
    eth_ret_freq = eth_resampled.pct_change().dropna()

    all_returns[f"{label} – BTC"] = btc_ret_freq
    all_returns[f"{label} – ETH"] = eth_ret_freq

diagnostics = dataFunc.test_distribution_diagnostics(all_returns)
print(diagnostics.round(3))


###########
# Univariate Models: In Sample
###########

returns = btc["Close"].resample("1D").last().pct_change().dropna()
results_table, best_model = funcUni.fit_garch_models(returns, dist="t")
print(results_table.sort_values("AIC").round(3))

param_table, info_table = funcUni.summarize_garch_model(best_model)
print("Model Parameters:")
print(param_table)
print("\nFit Statistics:")
print(info_table)

funcUni.plot_garch_residuals(best_model, title_prefix="BTC GARCH(1,3)")

# Autocorrelogram with corrected 95 ci???

# Threshold GARCH
# Regime Switching

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