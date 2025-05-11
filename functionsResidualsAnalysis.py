import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf

def robust_acf_plot(x, max_lag=None, title="Robust ACF"):
    x = np.asarray(x)
    n = len(x)
    if max_lag is None:
        max_lag = int(min(10 * np.log10(n), n - 1))

    def gamma(x, h):
        x = x - np.mean(x)
        h = abs(h)
        return np.sum(x[:n - h] * x[h:]) / n

    def rho(x, h):
        return gamma(x, h) / gamma(x, 0)

    acf_vals = [rho(x, h) for h in range(1, max_lag + 1)]
    x2 = x**2
    gamma0_sq = gamma(x, 0)**2
    variances = [(1 + gamma(x2, h) / gamma0_sq) / n for h in range(1, max_lag + 1)]
    bands = 1.96 * np.sqrt(variances)

    lower_band = -bands
    upper_band = bands

    minval = 1.2 * min(min(acf_vals), min(lower_band), -1.96 / np.sqrt(n))
    maxval = 1.2 * max(max(acf_vals), max(upper_band), 1.96 / np.sqrt(n))

    fig, ax = plt.subplots(figsize=(10, 4))
    lags = range(1, max_lag + 1)

    # Styled vertical lines
    for lag, val in zip(lags, acf_vals):
        ax.vlines(lag, 0, val, color='black', linewidth=1)

    # Dot markers
    ax.plot(lags, acf_vals, 'ko', markersize=4)

    # Corrected confidence bands
    ax.plot(lags, upper_band, linestyle='--', color='gray', linewidth=1)
    ax.plot(lags, lower_band, linestyle='--', color='gray', linewidth=1)
    ax.axhline(0, color='black', linewidth=0.5)

    ax.set_xlim(0.5, max_lag + 0.5)
    ax.set_ylim(minval, maxval)
    ax.set_xticks(np.arange(1, max_lag + 1))
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax.set_title(title)
    plt.tight_layout()
    plt.show(block=False)

