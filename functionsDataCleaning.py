import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

def compute_stats(r):
    return {
        "Average Return (%)": r.mean() * 100,
        "Annualized Volatility (%)": r.std() * np.sqrt(252 * (1440 / r.name)) * 100,
        "Skewness": r.skew(),
        "Excess Kurtosis": r.kurt()
    }

def load_close_series(path):
    df = pd.read_csv(path, parse_dates=["Open Time"], usecols=["Open Time", "Close"])
    df.set_index("Open Time", inplace=True)
    df["Close"] = pd.to_numeric(df["Close"], errors="coerce")
    return df

def plot_series(series_dict, title, ylabel):
    color_map = {
        "BTC": "blue",
        "ETH": "red"
    }
    plt.figure(figsize=(12, 4))
    for label, series in series_dict.items():
        plt.plot(series, label=label, color=color_map.get(label, None), alpha=0.6)
    plt.title(title)
    plt.xlabel("Date")
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show(block = False)


    # Intraday volatility → daily annualized
def daily_ann_vol(ret_1m):
    vol = ret_1m.pct_change().groupby(ret_1m.index.date).std()
    return pd.Series(vol.values * np.sqrt(1440 * 252) * 100, index=pd.to_datetime(vol.index))

def plot_autocorrelogram(series, lags=20, title="", squared=False):
    data = series ** 2 if squared else series
    T = len(data)
    conf_int = 1.96 / np.sqrt(T)

    fig, ax = plt.subplots(figsize=(10, 4))
    sm.graphics.tsa.plot_acf(data, lags=lags, alpha=None, ax=ax)
    ax.axhline(y=conf_int, linestyle='--', color='gray', linewidth=1)
    ax.axhline(y=-conf_int, linestyle='--', color='gray', linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("Lags")
    ax.set_ylabel("Autocorrelation")
    plt.tight_layout()
    plt.show(block=False)
