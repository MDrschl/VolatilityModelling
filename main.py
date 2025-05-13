import pandas as pd
import numpy as np
import os
import functionsCleaningDescriptives as dataFunc
import functionsUnivariate as funcUni
import matplotlib.pyplot as plt


from DCC_GARCH.GARCH.GARCH import GARCH
from DCC_GARCH.GARCH.GARCH_loss import garch_loss_gen
from DCC_GARCH.DCC.DCC import DCC
from DCC_GARCH.DCC.DCC_loss import dcc_loss_gen

###########
# Data Loading, Preprocessing, Descriptive Statistics
###########
btc_full = dataFunc.load_close_series("data/BTCUSDT_1m.csv").sort_index()
eth_full = dataFunc.load_close_series("data/ETHUSDT_1m.csv").sort_index()

split_date = "2023-10-31"

# Train set: up to and including 2023-10-31
btc_train = btc_full.loc[:split_date]
eth_train = eth_full.loc[:split_date]

returns_df = pd.concat([btc_ret, eth_ret], axis=1).dropna()
returns_df.columns = ["BTC", "ETH"]
returns_df.to_csv("train_returns.csv")

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

# Fit ARCH and GARCH models on 1-Day, 1-Hour and 6-Hour frequency

cryptocurrency = {
    "BTC": btc_train,
    "ETH": eth_train
}

frequency = {
    "1D": "Daily",
    "6H": "6-Hour",
    "1H": "Hourly"
    
}

for crypto, data in cryptocurrency.items():
    for freq, label in frequency.items():

        print(f"\n\n# {label} frequency analysis for {crypto}")
        print("-" * 50)

        # Extract return at specified frequency
        data_returns_freq = data["Close"].resample(freq).last().pct_change().dropna()

        # Fit ARCH and GARCH models and rank performance
        results, best_model = funcUni.fit_garch_models(data_returns_freq, dist="t")
        print(f"{crypto} GARCH Model Comparison ({label} frequency):\n", results.sort_values("AIC").round(3))

        # Return summary statistics on models
        params, info = funcUni.summarize_garch_model(best_model)
        print(f"\n{crypto} Model Parameters:\n", params)
        print(f"\n{crypto} Fit Statistics:\n", info)

        # Residual plots (Raw, Squared, and Standardized)
        funcUni.plot_garch_residuals(best_model, title_prefix=f"{crypto} Residuals ({label})")

        # Residual ACF (raw and squared)
        dataFunc.plot_autocorrelogram(best_model.resid.dropna(), title=f"{crypto} Residual ACF – {label}")
        dataFunc.plot_autocorrelogram(best_model.resid.dropna(), squared=True, title=f"{crypto} Squared Residual ACF – {label}")

        # Residual diagnostics table
        diagnostics = funcUni.garch_residual_diagnostics(best_model)
        print(f"\n{crypto} Residual Diagnostics:\n", diagnostics)






# Threshold GARCH

# Include PACF function


# Regime Switching

# Hourly, 6 hour frequency

###########
# Univariate Models: Out of Sample
###########

###########
# Multivariate Models: In Sample
###########

for freq, label in frequency.items():
    print(f"\n--- Topaceminem DCC-GARCH on {label} data ---")

    # Step 1: Resample and compute log returns
    btc = np.log(btc_train["Close"].resample(freq).last()).diff().dropna()
    eth = np.log(eth_train["Close"].resample(freq).last()).diff().dropna()

    # Step 2: Align both series
    combined = pd.concat([btc, eth], axis=1).dropna()
    combined.columns = ["BTC", "ETH"]

    # Step 3: Fit univariate GARCH for each asset
    btc_model = GARCH(1, 1)
    btc_model.set_loss(garch_loss_gen(1, 1))
    btc_model.set_max_itr(50)
    btc_model.fit(combined["BTC"])

    eth_model = GARCH(1, 1)
    eth_model.set_loss(garch_loss_gen(1, 1))
    eth_model.set_max_itr(50)
    eth_model.fit(combined["ETH"])

    # Step 4: Get standardized residuals
    btc_sigma = btc_model.sigma(combined["BTC"])
    eth_sigma = eth_model.sigma(combined["ETH"])
    epsilon = np.array([
        combined["BTC"] / btc_sigma,
        combined["ETH"] / eth_sigma
    ])

    # Step 5: Fit DCC model on standardized residuals
    model = DCC()
    model.set_loss(dcc_loss_gen())
    model.fit(epsilon)

    # Step 6: Plot dynamic correlation
    dynamic_corr = model.corr[:, 0, 1]
    plt.figure(figsize=(12, 5))
    plt.plot(combined.index[1:], dynamic_corr)  # Skip first obs due to lag
    plt.title(f'Dynamic Correlation BTC-ETH ({label})')
    plt.xlabel('Date')
    plt.ylabel('Correlation')
    plt.grid(True)
    plt.show()

# --- Topaceminem DCC-GARCH on Daily data ---

# Step 1: Resample and compute log returns
freq = '1D'
label = 'Daily'
btc = np.log(btc_train["Close"].resample(freq).last()).diff().dropna()
eth = np.log(eth_train["Close"].resample(freq).last()).diff().dropna()

# Step 2: Align both series
combined = pd.concat([btc, eth], axis=1).dropna()
combined.columns = ["BTC", "ETH"]

# Step 3: Fit univariate GARCH for each asset
btc_model = GARCH(1, 1)
btc_model.set_loss(garch_loss_gen(1, 1))
btc_model.set_max_itr(50)
btc_model.fit(combined["BTC"])
btc_model.get_theta()
eth_model = GARCH(1, 1)
eth_model.set_loss(garch_loss_gen(1, 1))
eth_model.set_max_itr(50)
eth_model.fit(combined["ETH"])
eth_model.get_theta()

# Step 4: Get standardized residuals
btc_sigma = btc_model.sigma(combined["BTC"])
eth_sigma = eth_model.sigma(combined["ETH"])
epsilon = np.array([
    combined["BTC"] / btc_sigma,
    combined["ETH"] / eth_sigma
])

# Step 5: Fit DCC model on standardized residuals
model = DCC()
model.set_loss(dcc_loss_gen())
model.fit(epsilon)

plt.figure(figsize=(12, 5))
plt.plot(combined.index[1:], dynamic_corr)  # Skip first obs due to lag
plt.title(f'Dynamic Correlation BTC-ETH ({label})')
plt.xlabel('Date')
plt.ylabel('Correlation')
plt.grid(True)
plt.show()



# VARMA-DCC-AGARCH
# Regime Switching

###########
# Multivariate Models: Out of Sample
###########


# Robust acf for codnitional mean, use robust
# Residuals univariate fitting 
#   Check flat bounds
# We already icnorporated hetehorskdeasityn when wemnestimaten garch model
