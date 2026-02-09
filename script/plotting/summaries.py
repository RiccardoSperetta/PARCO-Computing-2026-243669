# scripts/plotting/summaries.py
import pandas as pd
import numpy as np

def harmonic_mean(x):
    """Harmonic mean of positive values. Returns NaN if any <= 0."""
    if np.any(x <= 0):
        return np.nan
    return len(x) / np.sum(1.0 / x)


def compute_teps_harmonic(df: pd.DataFrame) -> pd.DataFrame:
    """
    harmonic mean per group (= same graph, method and number of cores)
    """
    return (
        df.groupby(['graph', 'variant', 'cores'])['TEPS']
          .agg(harmonic_mean)
          .reset_index(name='TEPS_harmonic')
    )


def compute_times_p90(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes 90th percentile for both total_time and comm_time (per group)
    """
    return (
        df.groupby(['graph', 'variant', 'cores'])
          .agg(
              total_time_p90=('total_time', lambda x: x.quantile(0.90)),
              comm_time_p90=('comm_time',   lambda x: x.quantile(0.90)),
          )
          .reset_index()
    )


def compute_max_over_mean_avg(df: pd.DataFrame) -> pd.DataFrame:
    """
    Mean of max_over_mean across runs (per group).
    """
    return (
        df.groupby(['graph', 'variant', 'cores'])['max_over_mean']
          .mean()  # or .median() if runs vary wildly
          .reset_index(name='max_over_mean_avg')
    )


def compute_cv_avg(df: pd.DataFrame) -> pd.DataFrame:
    """
    Mean of CV across runs (per group).
    """
    return (
        df.groupby(['graph', 'variant', 'cores'])['CV']
          .mean()  # or .median()
          .reset_index(name='CV_avg')
    )