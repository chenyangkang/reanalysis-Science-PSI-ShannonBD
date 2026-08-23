# Reanalysis of PSI and bird diversity

This repository contains two concise, self-contained reanalysis notebooks for the data from *China's solar expansion policy reduces bird diversity*.

- `01_replication_and_coefficient_reanalysis.ipynb` reconstructs the released analysis results presented in the paper, reproduces the principal fixed-effects and instrumental variable estimates, decomposes the reported R-squared values, checks the PSI-to-area pathway.
- `02_citizen_science_sampling_reanalysis.ipynb` examines the duration-per-observer variable, reporting discontinuities, sample inclusion, observer count, total observation time, sampling effort-adjusted outcomes, and IV sensitivity to observer count.

The notebooks have been executed from top to bottom and retain their outputs. Matching HTML exports are included for reading without Jupyter.

To see analysis and annotation, click on:
- [Notebook 1: Replication and coefficient reanalysis](01_replication_and_coefficient_reanalysis.ipynb)
- [Notebook 2: Citizen science sampling reanalysis](02_citizen_science_sampling_reanalysis.ipynb)

## Folder layout

Run the notebooks from this `scripts/` folder. They read the adjacent `../DATA/` folder and import only the local helper modules listed below:

- `load_data.py` builds the released analysis panel.
- `regression_utils.py` fits the fixed-effects models with PyFixest.
- `iv_observer_sensitivity.py` reproduces the IV model and observer-count sensitivity analysis.

The released data are read only; the notebooks do not modify files in `../DATA/`.

## Main findings

- The released negative PSI coefficient is reproduced, but PSI also predicts sampling efforts, including observer count and total recorded observation time within counties.
- The duration-per-observer variable does not preserve sampling scale, and the released effort fields contain extreme values, internal discrepancies, and strong sampling effort domain shift over time.
- Adding observer count changes the fixed-effects PSI estimate from negative and significant to small, positive, and nonsignificant. The corresponding richness estimate nearly disappears.
- The IV estimate is similarly eliminated by observer adjustment, while the instrument itself predicts observer count, raising concern about the exclusion restriction.
- These diagnostics do not prove that the ecological effect is zero. They show that the claimed negative causal effect cannot be distinguished from changing sampling coverage and detectability using the released aggregate data.

## Install

Python 3.11 was used for the saved results. From this folder, create an isolated environment and install the recorded dependencies with:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Run

Execute both notebooks and retain the new outputs:

```bash
jupyter nbconvert --to notebook --execute --inplace 01_replication_and_coefficient_reanalysis.ipynb
jupyter nbconvert --to notebook --execute --inplace 02_citizen_science_sampling_reanalysis.ipynb
```

Regenerate the readable HTML files after execution:

```bash
jupyter nbconvert --to html 01_replication_and_coefficient_reanalysis.ipynb
jupyter nbconvert --to html 02_citizen_science_sampling_reanalysis.ipynb
```

## Estimation notes

- The principal models absorb county and exact year-month fixed effects and cluster standard errors by county.
- Fixed-effects OLS uses PyFixest's recursive singleton removal. The authors' IV code specifies Stata's `keepsingleton`, so the IV reproduction uses `fixef_rm="none"`.
- The PyFixest small-sample correction is configured to match the released Stata specifications.

## DAGs
![DAGs](./figures/DAGs.png)

