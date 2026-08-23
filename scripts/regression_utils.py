from __future__ import annotations

import argparse
import base64
import html
import io
import math
import os
from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyfixest as pf
from scipy import stats

# Match the slope-only finite-sample correction used by the reported reghdfe reproductions.
PFIXEST_SSC = pf.ssc(k_adj=True, k_fixef="none", G_adj=True, G_df="min")

def _fe_formula(outcome: str, regressors: Sequence[str], fe_columns: Sequence[str]) -> str:
    formula = f"{outcome} ~ {' + '.join(regressors) if regressors else '1'}"
    if fe_columns:
        formula += f" | {' + '.join(fe_columns)}"
    return formula

def fit_fe(
    data: pd.DataFrame,
    outcome: str,
    regressors: Sequence[str],
    fe_columns: Sequence[str] = ("id", "year"), # although they named it year but this variable is actually year-month... The capitalization of the variable YEAR is the true year variable.
    cluster_column: str = "id",
    focal: str | None = "PI",
    label: str = "",
) -> dict:
    columns = list(dict.fromkeys([outcome, *regressors, *fe_columns, cluster_column]))
    sample = data[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    model = pf.feols(
        _fe_formula(outcome, regressors, fe_columns),
        data=sample,
        vcov={"CRV1": cluster_column},
        ssc=PFIXEST_SSC,
        fixef_rm="singleton",
        drop_intercept=not fe_columns,
    )

    coefficients = model.coef()
    standard_errors = model.se()
    residual = model.resid()
    ssr = float(residual @ residual)
    n = int(model._N)
    n_clusters = int(model._data[cluster_column].nunique())
    overall_r2 = float(model._r2)
    within_r2 = float(model._r2_within) if fe_columns else overall_r2
    fe_only_r2 = (
        1.0 - (1.0 - overall_r2) / (1.0 - within_r2)
        if fe_columns and within_r2 < 1.0
        else 0.0
    )

    partial_r2 = np.nan
    if focal is not None and focal in regressors and len(regressors) > 1:
        reduced_regressors = [name for name in regressors if name != focal] # leave this out
        reduced = pf.feols(
            _fe_formula(outcome, reduced_regressors, fe_columns),
            data=model._data,
            vcov={"CRV1": cluster_column},
            ssc=PFIXEST_SSC,
            fixef_rm="none",
            drop_intercept=not fe_columns,
        )
        reduced_residual = reduced.resid()
        reduced_ssr = float(reduced_residual @ reduced_residual)
        partial_r2 = (reduced_ssr - ssr) / reduced_ssr if reduced_ssr > 0 else np.nan

    result = {
        "label": label,
        "outcome": outcome,
        "fe": " + ".join(fe_columns),
        "n": n,
        "clusters": n_clusters,
        "overall_r2": overall_r2,
        "fe_only_r2": fe_only_r2,
        "within_r2": within_r2,
        "partial_r2_focal": partial_r2,
        "estimator": f"pyfixest {pf.__version__}",
        "ssr": ssr,
    }
    for name in regressors:
        coefficient = float(coefficients.get(name, np.nan))
        std_error = float(standard_errors.get(name, np.nan))
        p_value = 2 * stats.t.sf(abs(coefficient / std_error), df=max(1, n_clusters - 1))
        result[f"coef_{name}"] = coefficient
        result[f"se_{name}"] = std_error
        result[f"p_{name}"] = float(p_value)
    return result


