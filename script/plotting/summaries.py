# scripts/plotting/summaries.py
import pandas as pd
import numpy as np

def harmonic_mean(x):
    """Harmonic mean of positive values. Returns NaN if any <= 0."""
    if np.any(x <= 0):
        return np.nan
    return len(x) / np.sum(1.0 / x)


def compute_teps_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group by graph + variant + cores and compute:
    - harmonic mean of TEPS
    - arithmetic mean of total_time (for reference)
    - count of runs (should be ~64)
    """
    grouped = df.groupby(['graph', 'variant', 'cores'])

    summary = grouped.agg(
        TEPS_harmonic=('TEPS', harmonic_mean),
        TEPS_arithmetic=('TEPS', 'mean'),     # just for comparison / debugging
        total_time_mean=('total_time', 'mean'),
        count=('TEPS', 'count')
    ).reset_index()

    return summary