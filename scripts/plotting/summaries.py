# summaries.py
import pandas as pd


# ----------------------------------------------------------------------
# Helper:   Creates a summary with mean + 90th percentile for the given value columns.
#           Works for both time and perf DataFrames.
# ----------------------------------------------------------------------
def add_summary_stats(df: pd.DataFrame, value_cols: list):
    group_cols = ["matrix", "opt_level", "threads", "schedule", "chunksize"]
    
    agg_dict = {}
    for col in value_cols:
        agg_dict[col] = [
            ("mean", "mean"),
            ("p90",  lambda x: x.quantile(0.90))
        ]
    
    summary = df.groupby(group_cols, dropna=False).agg(agg_dict)
    # Flatten multi-level columns
    summary.columns = [f"{col}_{stat}" for col, stat in summary.columns]
    summary = summary.reset_index()
    
    return summary

# ----------------------------------------------------------------------
# Helper:   Adds GFLOPs and memory bandwidth...(?)
# ----------------------------------------------------------------------