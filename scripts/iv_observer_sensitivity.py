"""Reproduce Table S7 and test sensitivity to omitted observer count.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyfixest as pf
from scipy import stats


HERE = Path(__file__).resolve().parent
PROJECT_DIR = HERE.parent
sys.path.insert(0, str(HERE))

# import solar_bird_domain_shift_analysis as core  # noqa: E402
from load_data import build_analysis_data, winsorize
from regression_utils import fit_fe, PFIXEST_SSC

CONTROLS = [
    "Temp",
    "Wind",
    "Pop",
    "Duration",
    "Carbon",
    "Water",
    "Green",
    "Farm",
    "Grass",
]


def cluster_p_value(model, term: str, clusters: int) -> float:
    t_value = float(model.coef()[term] / model.se()[term])
    return float(2 * stats.t.sf(abs(t_value), df=max(1, clusters - 1)))


def build_iv_data(data_dir: Path) -> pd.DataFrame:
    df, _, _ = build_analysis_data(data_dir)
    sunshine = pd.read_stata(
        data_dir / "sunshine.dta", convert_categoricals=False
    )[["id", "YEAR", "sun"]].drop_duplicates(["id", "YEAR"])
    cpu = pd.read_stata(
        data_dir / "CCPU.dta", convert_categoricals=False
    )[["市", "YEAR", "CCPU"]].drop_duplicates(["市", "YEAR"])

    df = df.merge(sunshine, on=["id", "YEAR"], how="left", validate="many_to_one")
    df = df.merge(cpu, on=["市", "YEAR"], how="left", validate="many_to_one")

    # Match CODE.do: winsorize these source variables before forming the IV.
    df["sun"] = winsorize(df["sun"])
    df["CCPU"] = winsorize(df["CCPU"])
    df["Sun_ccpu"] = np.log(df["sun"]) / df["CCPU"]
    return df


def fit_iv(df: pd.DataFrame, add_observer_count: bool) -> dict[str, float | int | str]:
    extra = ["log_BN_w"] if add_observer_count else []
    columns = [
        "ShannonBD",
        "PI",
        "Sun_ccpu",
        *CONTROLS,
        *extra,
        "id",
        "year",
    ]
    sample = df[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()

    exogenous = [*CONTROLS, *extra]
    with warnings.catch_warnings():
        # PyFixest constructs an internal first-stage diagnostic that warns about
        # singletons even though the requested IV sample itself is retained.
        warnings.filterwarnings("ignore", message=r".*singleton fixed effect.*")
        iv_model = pf.feols(
            f"ShannonBD ~ {' + '.join(exogenous)} | id + year | PI ~ Sun_ccpu",
            data=sample,
            vcov={"CRV1": "id"},
            ssc=PFIXEST_SSC,
            fixef_rm="none",
        )
    first_stage = pf.feols(
        f"PI ~ Sun_ccpu + {' + '.join(exogenous)} | id + year",
        data=sample,
        vcov={"CRV1": "id"},
        ssc=PFIXEST_SSC,
        fixef_rm="none",
    )

    clusters = int(sample["id"].nunique())
    result: dict[str, float | int | str] = {
        "model": "IV + observer count" if add_observer_count else "Released Table S7 IV",
        "n": int(iv_model._N),
        "clusters": clusters,
        "estimator": f"pyfixest {pf.__version__}",
        "pi_coefficient": float(iv_model.coef()["PI"]),
        "pi_clustered_se": float(iv_model.se()["PI"]),
        "pi_p_value": cluster_p_value(iv_model, "PI", clusters),
        "first_stage_instrument_coefficient": float(first_stage.coef()["Sun_ccpu"]),
        "first_stage_instrument_clustered_se": float(first_stage.se()["Sun_ccpu"]),
        "first_stage_instrument_p_value": cluster_p_value(first_stage, "Sun_ccpu", clusters),
    }
    if add_observer_count:
        result.update(
            {
                "observer_coefficient": float(iv_model.coef()["log_BN_w"]),
                "observer_clustered_se": float(iv_model.se()["log_BN_w"]),
                "observer_p_value": cluster_p_value(iv_model, "log_BN_w", clusters),
            }
        )
    return result


def fit_instrument_observer_reduced_form(df: pd.DataFrame) -> dict[str, float | int | str]:
    columns = ["log_BN_w", "Sun_ccpu", *CONTROLS, "id", "year"]
    sample = df[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    model = pf.feols(
        f"log_BN_w ~ Sun_ccpu + {' + '.join(CONTROLS)} | id + year",
        data=sample,
        vcov={"CRV1": "id"},
        ssc=PFIXEST_SSC,
        fixef_rm="none",
    )
    clusters = int(sample["id"].nunique())
    return {
        "model": "Instrument predicts observer count",
        "n": int(model._N),
        "clusters": clusters,
        "estimator": f"pyfixest {pf.__version__}",
        "pi_coefficient": np.nan,
        "pi_clustered_se": np.nan,
        "pi_p_value": np.nan,
        "first_stage_instrument_coefficient": float(model.coef()["Sun_ccpu"]),
        "first_stage_instrument_clustered_se": float(model.se()["Sun_ccpu"]),
        "first_stage_instrument_p_value": cluster_p_value(model, "Sun_ccpu", clusters),
    }
