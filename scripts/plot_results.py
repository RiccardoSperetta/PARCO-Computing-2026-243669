import os
from pathlib import Path
from typing import List, Dict
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(".")                     # project root
RESULTS_DIR = PROJECT_ROOT / "results"

# Expected values (adjust only if your actual folders differ)
OPT_LEVELS = ["O0", "O1", "O2", "O3", "Ofast"]
THREAD_COUNTS = [2, 4, 8, 16, 32, 64, 128]
SCHEDULES = ["static", "dynamic", "guided"]
CHUNKSIZES = [10, 100, 1000]

# ----------------------------------------------------------------------
# Helper: read time.txt (values are in milliseconds, in scientific notation)
# ----------------------------------------------------------------------
def read_times_ms(time_file: Path) -> List[float]:
    """Read time.txt and return list of floats in milliseconds."""
    if not time_file.exists():
        return []
    
    # Read as string → handles 1.23e-02 perfectly
    raw = pd.read_csv(time_file, header=None, dtype=str).iloc[:, 0].str.strip()
    times_ms = pd.to_numeric(raw, errors='coerce').dropna()
    return times_ms.tolist()

# ----------------------------------------------------------------------
# Helper: read perf.txt (10 lines of "L1_miss% LLC_miss%")
# ----------------------------------------------------------------------
def read_perf(perf_file: Path) -> List[Dict]:
    """Return list of dicts with L1 and LLC miss percentages."""
    records = []
    if not perf_file.exists():
        return records
    
    with open(perf_file) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            l1, llc = map(float, parts)
            records.append({"run": i, "L1_miss_percent": l1, "LLC_miss_percent": llc})
    return records


# ----------------------------------------------------------------------
# Helper: read perf.txt (10 lines of "L1_miss% LLC_miss%")
# ----------------------------------------------------------------------
def add_summary_stats(df: pd.DataFrame, value_cols: list):
    """
    Creates a summary with mean + 90th percentile for the given value columns.
    Works for both time and perf DataFrames.
    """
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



### --------------------------------------------------------------------
# "MAIN"
### --------------------------------------------------------------------

# ----------------------------------------------------------------------
# Storage for all raw measurements
# ----------------------------------------------------------------------
time_records: List[Dict] = []
perf_records: List[Dict] = []

# ----------------------------------------------------------------------
# Walk through all matrices
# ----------------------------------------------------------------------
for matrix_path in RESULTS_DIR.iterdir():
    if not matrix_path.is_dir():
        continue
    matrix = matrix_path.name
    print(f"Processing matrix: {matrix}")

    # ====================== SEQUENTIAL RUNS ======================
    seq_path = matrix_path / "sequential"
    if seq_path.exists():
        for opt_folder in seq_path.iterdir():
            if not opt_folder.is_dir() or opt_folder.name not in OPT_LEVELS:
                continue
            opt = opt_folder.name

            time_file = opt_folder / "time.txt"
            perf_file = opt_folder / "perf.txt"

            # --- times ---
            for i, t_ms in enumerate(read_times_ms(time_file), 1):
                time_records.append({
                    "matrix": matrix,
                    "opt_level": opt,
                    "threads": 1,
                    "schedule": "sequential",
                    "chunksize": None,
                    "run": i,
                    "time_ms": t_ms                  
                })

            # --- perf ---
            for rec in read_perf(perf_file):
                perf_records.append({
                    "matrix": matrix,
                    "opt_level": opt,
                    "threads": 1,
                    "schedule": "sequential",
                    "chunksize": None,
                    **rec # includes run, L1 and LLC cache misses %
                })

    # ====================== PARALLEL RUNS ======================
    for th_folder in matrix_path.iterdir():
        if not th_folder.is_dir() or not th_folder.name.startswith("th"):
            continue
        try:
            threads = int(th_folder.name[2:])
            if threads not in THREAD_COUNTS:
                continue
        except ValueError:
            continue

        for config_folder in th_folder.iterdir():
            if not config_folder.is_dir():
                continue
            folder_name = config_folder.name
            if "_" not in folder_name:
                continue

            schedule, chunksize_str = folder_name.rsplit("_", 1)
            if schedule not in SCHEDULES:
                continue
            try:
                chunksize = int(chunksize_str)
                if chunksize not in CHUNKSIZES:
                    continue
            except ValueError:
                continue

            time_file = config_folder / "time.txt"
            perf_file = config_folder / "perf.txt"

            # --- times ---
            for i, t_ms in enumerate(read_times_ms(time_file), 1):
                time_records.append({
                    "matrix": matrix,
                    "opt_level": "O3",               
                    "threads": threads,
                    "schedule": schedule,
                    "chunksize": chunksize,
                    "run": i,
                    "time_ms": t_ms                    
                })

            # --- perf ---
            for rec in read_perf(perf_file):
                perf_records.append({
                    "matrix": matrix,
                    "opt_level": "O3",
                    "threads": threads,
                    "schedule": schedule,
                    "chunksize": chunksize,
                    **rec   # includes run, L1 and LLC cache misses %
                })

# ----------------------------------------------------------------------
# Build final DataFrames
# ----------------------------------------------------------------------
df_time = pd.DataFrame(time_records)
df_perf = pd.DataFrame(perf_records)


# ====================== SUMMARIES WITH 90TH PERCENTILE ======================
df_time_summary = add_summary_stats(df_time, ["time_ms"])
df_perf_summary = add_summary_stats(df_perf, ["L1_miss_percent", "LLC_miss_percent"])

# debug prints
print("\n--- df_time_summary (first 20 rows) ---")
print(df_time_summary.head(20).to_string(index=False))

print("\n--- df_perf_summary (first 15 rows) ---")
print(df_perf_summary.head(15).to_string(index=False))



# ====================== PLOTTING ======================

for matrix_path in RESULTS_DIR.iterdir():
    if not matrix_path.is_dir():
        continue
    matrix = matrix_path.name
    SCHEDULE = "static"
    CHUNKSIZE = 100

    # Create output folder
    plot_dir = Path("plots") / matrix
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Filter data: FIRST sequential VS parallel (stati, 100)
    seq = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["threads"] == 1)
    ].copy()

    par = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["schedule"] == SCHEDULE) &
        (df_time_summary["chunksize"] == CHUNKSIZE)
    ].copy()

    # Build nice labels and correct order
    seq["label"] = seq["opt_level"]
    seq["order"] = seq["opt_level"].map({"O0":0, "O1":1, "O2":2, "O3":3, "Ofast":4})

    par["label"] = par["threads"].astype(str) + " th"
    par["order"] = par["threads"].map({2:5, 4:6, 8:7, 16:8, 32:9, 64:10, 128:11})

    # sequential and parallel code confronted:
    data = pd.concat([seq, par], ignore_index=True)
    data = data.sort_values("order")

    # Plot
    if data.empty:
        print(f"No data to plot for {matrix} ({SCHEDULE} chunksize {CHUNKSIZE})")
    else:
        plt.figure(figsize=(11, 6))

        # ensure non-negative yerr (in case p90 < mean due to small samples) and convert to numpy
        yerr = (data["time_ms_p90"] - data["time_ms_mean"]).clip(lower=0).values

        # build color list and ensure it matches the length of the data
        colors = ["#ff9999"] * len(seq) + ["#66b3ff"] * len(par)
        if len(colors) < len(data):
            colors = (colors + ["#66b3ff"] * len(data))[:len(data)]
        else:
            colors = colors[:len(data)]

        bars = plt.bar(
            data["label"],
            data["time_ms_mean"],
            yerr=yerr,
            capsize=6,
            color=colors,
            edgecolor="black",
            linewidth=0.8
        )

        plt.title(f"{matrix} — {SCHEDULE}_chunksize{CHUNKSIZE}\nMean time with 90th percentile error bars", fontsize=14)
        plt.ylabel("Execution time (ms)", fontsize=12)
        plt.xlabel("Configuration", fontsize=12)
        #plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        # Save automatically
        filename = plot_dir / f"time_{SCHEDULE}_{CHUNKSIZE}.png"
        plt.savefig(filename, dpi=300, bbox_inches="tight")

        plt.close()

        print(f"Plot saved → {filename}")
