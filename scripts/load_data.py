from __future__ import annotations

import argparse
import base64
import html
import io
import math
import os
from pathlib import Path
from typing import Sequence
import numpy as np
import pandas as pd

def read_dta(data_dir: Path, name: str) -> pd.DataFrame:
    return pd.read_stata(data_dir / name, convert_categoricals=False)



def winsorize(s: pd.Series, lower: float = 0.01, upper: float = 0.99) -> pd.Series:
    nonmissing = s.dropna()
    if nonmissing.empty:
        return s
    lo, hi = nonmissing.quantile([lower, upper]).to_numpy()
    return s.clip(lo, hi)



def checked_left_merge(
    master: pd.DataFrame,
    using: pd.DataFrame,
    keys: Sequence[str],
    columns: Sequence[str],
    validate: str,
) -> pd.DataFrame:
    right = using[list(keys) + list(columns)].copy()
    return master.merge(right, how="left", on=list(keys), validate=validate)




def build_analysis_data(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    bird = read_dta(data_dir, "Bird.dta")
    bird["YEAR"] = bird["YEAR"].astype(int)
    bird["month"] = bird["month"].astype(int)
    bird["year"] = pd.to_datetime(bird["year"])

    btn = read_dta(data_dir, "BTN.dta")
    btn["year"] = pd.to_datetime(btn["year"])
    df = checked_left_merge(bird, btn, ["id", "year"], ["BT", "BN"], "one_to_one")

    annual_sources = [
        ("PI.dta", ["PI"], {}),
        ("Climate.dta", ["Avt", "Wind"], {"Avt": "Temp"}),
        ("PD.dta", ["PD"], {}),
        ("air.dta", ["CO2", "PM"], {"PM": "PM25"}),
        ("TD.dta", ["Water", "Green", "Farm", "Grass"], {}),
        ("IC.dta", ["ICtotal"], {}),
        ("Ndvi.dta", ["Ndvi"], {}),
        ("NL.dta", ["Nightlight"], {}),
        ("LAI.dta", ["LAI"], {}),
    ]
    for filename, columns, rename in annual_sources:
        source = read_dta(data_dir, filename).rename(columns=rename)
        use_columns = [rename.get(c, c) for c in columns]
        if "YEAR" in source:
            source["YEAR"] = source["YEAR"].astype(int)
        df = checked_left_merge(df, source, ["id", "YEAR"], use_columns, "many_to_one")

    richness = read_dta(data_dir, "Richeven.dta")
    richness["year"] = pd.to_datetime(richness["year"])
    df = checked_left_merge(df, richness, ["id", "year"], ["丰富度"], "one_to_one")
    df = df.rename(columns={"丰富度": "Richness"})

    # Reproduce the released Stata transformations.
    df["PI_missing_in_source"] = df["PI"].isna().astype(int)
    df["PI"] = df["PI"].fillna(0.0)
    df["ICtotal"] = df["ICtotal"].fillna(0.0)
    df["BTN"] = df["BT"] / df["BN"]

    required = ["Temp", "Wind", "PD", "BTN", "CO2", "Water", "Green", "Farm", "Grass"]
    df = df.dropna(subset=required).copy()
    df = df[(df["PD"] > 0) & (df["CO2"] > 0) & (df["BN"] > 0)].copy()

    df["Area"] = df["ICtotal"] * 1000.0 / 0.15 / 1_000_000.0
    df["Pop"] = np.log(df["PD"])
    df["Duration"] = np.log1p(df["BTN"])
    df["Carbon"] = np.log(df["CO2"])
    df["PM"] = np.log(df["PM25"])
    df["log_BN"] = np.log1p(df["BN"])
    df["log_BT"] = np.log1p(df["BT"])
    df["BT_zero"] = (df["BT"] == 0).astype(float)
    df["Ndvi"] = df["Ndvi"] * 100.0
    df["Evenness"] = np.where(df["Richness"] > 1, df["ShannonBD"] / np.log(df["Richness"]), np.nan)
    df["county_month"] = df["id"].astype(str) + "_" + df["month"].astype(str).str.zfill(2)
    df["trend"] = (df["year"].dt.year - 2014) * 12 + df["year"].dt.month - 1

    released_winsor_vars = [
        "ShannonBD", "SimpsonBD", "PI", "Temp", "Wind", "Pop", "Duration", "Carbon",
        "Water", "Green", "Farm", "Grass", "PM", "ICtotal", "Area", "Ndvi",
        "Nightlight", "LAI", "Richness", "Evenness",
    ]
    for col in released_winsor_vars:
        if col in df:
            df[col] = winsorize(df[col])
    # Diagnostics not used in the paper: control extreme effort values without altering source columns.
    df["log_BN_w"] = winsorize(df["log_BN"])
    df["log_BT_w"] = winsorize(df["log_BT"])

    # Policy panel for constructing leads and the full county-month observation grid.
    pi = read_dta(data_dir, "PI.dta")[["id", "YEAR", "PI"]].copy()
    pi["YEAR"] = pi["YEAR"].astype(int)
    pi = pi.sort_values(["id", "YEAR"])
    for k in (1, 2):
        pi[f"PI_lead{k}"] = pi.groupby("id")["PI"].shift(-k)
        pi[f"PI_lag{k}"] = pi.groupby("id")["PI"].shift(k)
    df = checked_left_merge(df, pi, ["id", "YEAR"], ["PI_lead1", "PI_lead2", "PI_lag1", "PI_lag2"], "many_to_one")
    for col in ["PI_lead1", "PI_lead2", "PI_lag1", "PI_lag2"]:
        df[col] = winsorize(df[col])

    return df.sort_values(["id", "year"]).reset_index(drop=True), bird, pi


def make_full_observation_grid(
    analysis: pd.DataFrame,
    bird: pd.DataFrame,
    pi: pd.DataFrame,
) -> pd.DataFrame:
    ids = np.sort(analysis["id"].unique())
    dates = pd.date_range("2014-01-01", "2023-12-01", freq="MS")
    grid = pd.MultiIndex.from_product([ids, dates], names=["id", "year"]).to_frame(index=False)
    grid["YEAR"] = grid["year"].dt.year.astype(int)
    observed = bird[["id", "year"]].drop_duplicates().assign(observed=1)
    observed["year"] = pd.to_datetime(observed["year"])
    grid = grid.merge(observed, how="left", on=["id", "year"], validate="one_to_one")
    grid["observed"] = grid["observed"].fillna(0.0)
    grid = grid.merge(pi[["id", "YEAR", "PI"]], how="left", on=["id", "YEAR"], validate="many_to_one")
    grid["PI"] = grid["PI"].fillna(0.0)
    grid["PI"] = winsorize(grid["PI"])
    return grid