import pandas as pd
import numpy as np
from arch import arch_model
import matplotlib.pyplot as plt

def fit_garch_models(returns, model_name="Series", dist="normal", max_order=3):
    """
    Fits ARCH(p) and GARCH(p,q) models for p, q in 1..max_order and finds best model by AIC.
    
    Parameters:
    - returns: pd.Series of returns
    - model_name: str label for output
    - dist: error distribution ('normal', 't', etc.)
    - max_order: max lag order for p and q (default: 3)
    
    Returns:
    - pd.DataFrame with model diagnostics (LLH, AIC, BIC, params)
    - best_model: model object with lowest AIC
    """
    results = []
    best_model = None
    best_aic = np.inf  # ← initialize this

    # ARCH(p)
    for p in range(1, max_order + 1):
        try:
            model = arch_model(returns, vol='GARCH', p=p, q=0, dist=dist, rescale=True)
            res = model.fit(disp="off")
            results.append({
                "Model": f"ARCH({p})",
                "Log-Likelihood": res.loglikelihood,
                "AIC": res.aic,
                "BIC": res.bic,
                "Params": res.params.to_dict()
            })
            if res.aic < best_aic:
                best_aic = res.aic
                best_model = res
        except Exception as e:
            results.append({
                "Model": f"ARCH({p})",
                "Log-Likelihood": None,
                "AIC": None,
                "BIC": None,
                "Params": str(e)
            })

    # GARCH(p,q)
    for p in range(1, max_order + 1):
        for q in range(1, max_order + 1):
            try:
                model = arch_model(returns, vol='GARCH', p=p, q=q, dist=dist, rescale=True)
                res = model.fit(disp="off")
                results.append({
                    "Model": f"GARCH({p},{q})",
                    "Log-Likelihood": res.loglikelihood,
                    "AIC": res.aic,
                    "BIC": res.bic,
                    "Params": res.params.to_dict()
                })
                if res.aic < best_aic:
                    best_aic = res.aic
                    best_model = res
            except Exception as e:
                results.append({
                    "Model": f"GARCH({p},{q})",
                    "Log-Likelihood": None,
                    "AIC": None,
                    "BIC": None,
                    "Params": str(e)
                })

    return pd.DataFrame(results), best_model



def plot_garch_residuals(fitted_model, lags=20, title_prefix=""):
    """
    Plots residuals, squared residuals, and standardized residuals from a fitted GARCH model.
    
    Parameters:
    - fitted_model: result from arch_model(...).fit()
    - lags: number of lags for autocorrelation plots
    - title_prefix: optional prefix for plot titles
    """
    resid = fitted_model.resid
    std_resid = fitted_model.std_resid
    sq_resid = resid ** 2

    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axs[0].plot(resid)
    axs[0].set_title(f"{title_prefix} Raw Residuals")

    axs[1].plot(sq_resid)
    axs[1].set_title(f"{title_prefix} Squared Residuals")

    axs[2].plot(std_resid)
    axs[2].set_title(f"{title_prefix} Standardized Residuals")

    for ax in axs:
        ax.grid(True)

    plt.tight_layout()
    plt.show(block = False)


def summarize_garch_model(fitted_model):
    """
    Returns two tables:
    1. Parameter estimates with standard errors, t-stats, and p-values.
    2. Model fit statistics (Log-Likelihood, AIC, BIC).

    Parameters:
    - fitted_model: arch.univariate.base.ARCHModelResult

    Returns:
    - param_table: pd.DataFrame with parameter diagnostics
    - info_table: pd.DataFrame with model-level stats
    """
    params = fitted_model.params
    stderr = fitted_model.std_err
    tstat = params / stderr
    pval = fitted_model.pvalues

    param_table = pd.DataFrame({
        "Parameter": params,
        "Std. Error": stderr,
        "t-Statistic": tstat,
        "p-Value": pval
    })

    info = pd.DataFrame({
        "Log-Likelihood": [fitted_model.loglikelihood],
        "AIC": [fitted_model.aic],
        "BIC": [fitted_model.bic]
    })

    return param_table.round(4), info.round(4)


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
    robust_bands = 1.96 * np.sqrt(variances)

    # Standard (non-robust) CI
    standard_ci = 1.96 / np.sqrt(n)

    lower_robust = -robust_bands
    upper_robust = robust_bands

    minval = 1.2 * min(min(acf_vals), min(lower_robust), -standard_ci)
    maxval = 1.2 * max(max(acf_vals), max(upper_robust), standard_ci)

    fig, ax = plt.subplots(figsize=(10, 4))
    lags = range(1, max_lag + 1)

    # Vertical lines and dots
    for lag, val in zip(lags, acf_vals):
        ax.vlines(lag, 0, val, color='black', linewidth=1)
    ax.plot(lags, acf_vals, 'ko', markersize=4)

    # Robust CI
    ax.plot(lags, upper_robust, linestyle='--', color='gray', linewidth=1)
    ax.plot(lags, lower_robust, linestyle='--', color='gray', linewidth=1)

    # Standard CI
    ax.axhline(standard_ci, linestyle='--', color='red', linewidth=1, label='Standard 95% CI')
    ax.axhline(-standard_ci, linestyle='--', color='red', linewidth=1)


    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xlim(0.5, max_lag + 0.5)
    ax.set_ylim(minval, maxval)
    ax.set_xticks(np.arange(1, max_lag + 1))
    ax.set_xlabel("Lag")
    ax.set_ylabel("Autocorrelation")
    ax.set_title(title)
    ax.plot([], [], 'k--', label='Robust 95% CI')
    ax.plot([], [], 'r:', label='Standard 95% CI')
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.show(block=False)
