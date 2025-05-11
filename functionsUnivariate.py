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
