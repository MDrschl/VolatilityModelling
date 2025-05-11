import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf


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
    acf_vals = acf(data, nlags=lags, fft=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    markerline, stemlines, baseline = ax.stem(range(len(acf_vals)), acf_vals)
    plt.setp(markerline, markersize=4, marker='o', color='black')
    plt.setp(stemlines, linewidth=1.2, color='black')
    plt.setp(baseline, linewidth=0)

    ax.axhline(y=0, linestyle='-', color='black', linewidth=0.5)
    ax.axhline(y=conf_int, linestyle='--', color='gray', linewidth=1)
    ax.axhline(y=-conf_int, linestyle='--', color='gray', linewidth=1)

    ax.set_title(title)
    ax.set_xlabel("Lags")
    ax.set_ylabel("Autocorrelation")
    plt.tight_layout()
    plt.show(block=False)

import pandas as pd
import numpy as np
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import jarque_bera

def test_distribution_diagnostics(series_dict, lags=20):
    """
    Performs Ljung-Box test for autocorrelation, McLeod-Li test on squared returns, 
    and Jarque-Bera test for normality across return series.

    Parameters:
    - series_dict: dict with {label: return_series}
    - lags: number of lags for Ljung-Box/McLeod-Li (default 20)

    Returns:
    - pd.DataFrame with test stats and p-values for each series
    """
    results = []

    for label, series in series_dict.items():
        series = series.dropna()
        squared = series ** 2

        # Ljung-Box
        lb = acorr_ljungbox(series, lags=[lags], return_df=True)
        lb_stat, lb_p = lb.iloc[0]['lb_stat'], lb.iloc[0]['lb_pvalue']

        # McLeod-Li
        ml = acorr_ljungbox(squared, lags=[lags], return_df=True)
        ml_stat, ml_p = ml.iloc[0]['lb_stat'], ml.iloc[0]['lb_pvalue']

        # Jarque-Bera
        jb_stat, jb_p, _, _ = jarque_bera(series)

        results.append({
            "Series": label,
            "Ljung-Box stat": lb_stat,
            "Ljung-Box p": lb_p,
            "McLeod-Li stat": ml_stat,
            "McLeod-Li p": ml_p,
            "Jarque-Bera stat": jb_stat,
            "Jarque-Bera p": jb_p
        })

    return pd.DataFrame(results)
