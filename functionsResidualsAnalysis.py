import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf


def garch_corrected_acf_ci(residuals, lags=20, alpha=0.05):
    """
    Computes GARCH-corrected confidence intervals for sample autocorrelations.
    
    Parameters:
    - residuals: array-like, residuals or demeaned returns
    - lags: number of lags to include
    - alpha: significance level (default = 0.05 for 95% CI)

    Returns:
    - acf_vals: np.ndarray of autocorrelations
    - ci_bounds: np.ndarray of symmetric confidence bounds per lag
    """
    residuals = np.asarray(residuals)
    T = len(residuals)
    eps = residuals - np.mean(residuals)
    gamma0 = np.mean(eps ** 2)

    acf_vals = acf(eps, nlags=lags, fft=True)
    var_rho = []

    for j in range(lags + 1):
        if j == 0:
            var_rho.append(0.0)
            continue
        terms = eps[j:] ** 2 * eps[:-j] ** 2
        sigma_gamma_jj = np.mean(terms)
        var_rho_j = sigma_gamma_jj / (gamma0 ** 2) / T
        var_rho.append(var_rho_j)

    z = 1.96  # 95% confidence level
    ci_bounds = z * np.sqrt(var_rho)

    return acf_vals, ci_bounds

def plot_garch_corrected_acf(series, lags=20, title="GARCH-Corrected Autocorrelogram"):
    acf_vals, ci_bounds = garch_corrected_acf_ci(series, lags=lags)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.stem(range(len(acf_vals)), acf_vals, markerfmt='ko', basefmt=" ", linefmt='k-', bottom=0)
    ax.axhline(y=0, linestyle='-', color='black', linewidth=0.5)

    for lag in range(1, len(acf_vals)):
        ax.axhline(y=ci_bounds[lag], linestyle='--', color='gray', linewidth=1)
        ax.axhline(y=-ci_bounds[lag], linestyle='--', color='gray', linewidth=1)

    ax.set_title(title)
    ax.set_xlabel("Lags")
    ax.set_ylabel("Autocorrelation")
    plt.tight_layout()
    plt.show(block = False)
