import pandas as pd
import numpy as np
import functionsDataCleaning as dataFunc

btc = dataFunc.load_close_series("data/BTCUSDT_1m.csv")
eth = dataFunc.load_close_series("data/ETHUSDT_1m.csv")

# Resample to daily prices and returns
btc_daily = btc["Close"].resample("1D").last().dropna()
eth_daily = eth["Close"].resample("1D").last().dropna()
btc_ret = btc_daily.pct_change().dropna()
eth_ret = eth_daily.pct_change().dropna()

# Plot 1: Daily prices
dataFunc.plot_series({"BTC": btc_daily, "ETH": eth_daily}, "Daily Close Price", "Price (USDT)")

# Plot 2: Daily returns
dataFunc.plot_series({"BTC": btc_ret, "ETH": eth_ret}, "Daily Returns", "Return")

# Summary stats
def compute_stats(r):
    return {
        "Average Return (%)": r.mean() * 100,
        "Annualized Volatility (%)": r.std() * np.sqrt(252) * 100,
        "Skewness": r.skew(),
        "Excess Kurtosis": r.kurt()
    }

summary = pd.DataFrame([compute_stats(btc_ret), compute_stats(eth_ret)], index=["BTC", "ETH"])
print(summary.round(3))


vol_btc = dataFunc.daily_ann_vol(btc["Close"].dropna())
vol_eth = dataFunc.daily_ann_vol(eth["Close"].dropna())

# Plot 3: Daily annualized volatility
dataFunc.plot_series({"BTC": vol_btc, "ETH": vol_eth}, "Daily Annualized Volatility (from 1-Minute Data)", "Volatility (%)")
