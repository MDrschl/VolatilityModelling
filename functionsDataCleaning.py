import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and preprocess
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
        plt.plot(series, label=label, color=color_map.get(label, None), alpha=0.3)
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
