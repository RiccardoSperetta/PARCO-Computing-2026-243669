import os
from pathlib import Path
from typing import List, Dict
import pandas as pd
import matplotlib.pyplot as plt

from loaders import load_all_data
from summaries import add_summary_stats



### --------------------------------------------------------------------
# PLOTTING FUNCTIONS
### --------------------------------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def plot_time_line(matrix: str, df_time_summary: pd.DataFrame, plot_dir: Path,
                   schedule: str = "static", chunksize: int = 100):
    """Line graph of execution time (mean) with 90th percentile as the actual points"""
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Sequential O0–Ofast
    seq = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["threads"] == 1)
    ].copy()
    seq["label"] = seq["opt_level"]
    seq["order"] = seq["opt_level"].map({"O0":0, "O1":1, "O2":2, "O3":3, "Ofast":4})

    # Parallel
    par = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["schedule"] == schedule) &
        (df_time_summary["chunksize"] == chunksize)
    ].copy()
    par["label"] = par["threads"].astype(str) + " th"
    par["order"] = par["threads"].map({2:5,4:6,8:7,16:8,32:9,64:10,128:11})


    data = pd.concat([seq, par]).sort_values("order")

    plt.figure(figsize=(10, 5.5))
    # Line + points follow the 90th percentile
    plt.plot(data["label"], data["time_ms_p90"], 
             marker='o', markersize=7, linewidth=2.5, color="#1f77b4")
    plt.fill_between(data["label"], 
                     data["time_ms_mean"], 
                     data["time_ms_p90"], 
                     alpha=0.15, color="#1f77b4")

    plt.title(f"Execution time — {matrix}\n{schedule}, chunksize = {chunksize}", fontsize=14)
    plt.ylabel("Time (ms) — 90th percentile", fontsize=12)
    plt.xlabel("Configuration")
    plt.xticks(rotation=45)
    plt.ylim(0, None)
    plt.tight_layout()

    plt.savefig(plot_dir / f"time_line_{schedule}_{chunksize}.png", dpi=300)
    plt.close()


def plot_cache_misses_bar(matrix: str, df_perf_summary: pd.DataFrame, plot_dir: Path):
    """Grouped bar chart: L1 and LLC miss % for O3 sequential + parallel threads (static_100 only"""
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Keep only O3 sequential + static_100 parallel
    data = df_perf_summary[
        (df_perf_summary["matrix"] == matrix) &
        (
            ((df_perf_summary["threads"] == 1) & (df_perf_summary["opt_level"] == "O3")) |
             (df_perf_summary["schedule"] == "static") & (df_perf_summary["chunksize"] == 100)
        )
    ].copy()

    if data.empty:
        return

    # assigns correct name to each set of bars:
    data["config"] = data.apply(
        lambda row: "O3 (seq)" if row["threads"] == 1 else f"{int(row['threads'])} th",
        axis=1
    )
    data = data.sort_values("threads")

    x = range(len(data))
    width = 0.35

    # shows both L1 and LLc cache misses % (90% percentile) one next to the other
    plt.figure(figsize=(11, 5.5))
    plt.bar([i - width/2 for i in x], data["L1_miss_percent_p90"],
            width, label="L1 miss %", color="#ff7f0e", edgecolor="black")
    plt.bar([i + width/2 for i in x], data["LLC_miss_percent_p90"],
            width, label="LLC miss %", color="#2ca02c", edgecolor="black")

    plt.xticks(x, data["config"], rotation=45)
    plt.ylabel("Cache miss rate — 90th percentile (%)")
    plt.title(f"Cache miss rates — {matrix}\n(static schedule, chunksize=100)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(plot_dir / "cache_misses.png", dpi=300)
    plt.close()


def plot_scheduling_comparison(matrix: str, df_time_summary: pd.DataFrame, plot_dir: Path):
    threads_list = [8, 16, 32, 64, 128]
    """
    Suggested layout (you can choose one):
    Option A (my favourite) → 3 × 3 subplots (schedule × chunksize)
    Option B → one big grouped bar per thread count (15 groups per thread)
    """
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Keep only parallel runs with ≥8 threads
    data = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["threads"] >= 8)
    ].copy()

    if data.empty:
        return

    threads = sorted(data["threads"].unique())
    schedules = ["static", "dynamic", "guided"]
    chunksizes = [10, 100, 1000]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    bar_width = 0.25
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    for ax, chunksize in zip(axes, chunksizes):
        subset = data[data["chunksize"] == chunksize]

        for i, sched in enumerate(schedules):
            sched_data = subset[subset["schedule"] == sched]
            # create bars at position = thread_index + small offset
            x_pos = [threads.index(t) + i*bar_width for t in sched_data["threads"]]
            ax.bar(x_pos,
                   sched_data["time_ms_mean"],
                   width=bar_width,
                   label=sched if ax == axes[0] else "",
                   color=colors[i],
                   edgecolor="black")

        ax.set_title(f"chunksize = {chunksize}")
        ax.set_xlabel("Threads")
        ax.set_xticks([threads.index(t) + bar_width for t in threads])
        ax.set_xticklabels(threads)

    axes[0].set_ylabel("Time (ms) — arithmetic mean")
    fig.suptitle(f"Scheduling comparison — {matrix}", fontsize=16)
    fig.legend(schedules, loc="upper center", ncol=3)
    plt.tight_layout()
    plt.savefig(plot_dir / "scheduling_comparison.png", dpi=300)
    plt.close()


### --------------------------------------------------------------------
# "MAIN"
### --------------------------------------------------------------------
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
# Helper:   From the loaders.py "load_all_data" function I retrieve data frames loaded 
#           with all the info stored in the results folder
# ----------------------------------------------------------------------
[df_time, df_perf] = load_all_data(RESULTS_DIR)


# ====================== SUMMARIES WITH 90TH PERCENTILE ======================
# now the previous data frames are summarised, showing their mean and 90% percentile
df_time_summary = add_summary_stats(df_time, ["time_ms"])
df_perf_summary = add_summary_stats(df_perf, ["L1_miss_percent", "LLC_miss_percent"])


# debug prints
print("\n--- df_time_summary (first 20 rows) ---")
print(df_time_summary.head(20).to_string(index=False))

print("\n--- df_perf_summary (first 15 rows) ---")
print(df_perf_summary.head(15).to_string(index=False))


# iterative call on all matrices to build plots:
for matrix_path in RESULTS_DIR.iterdir():
    if not matrix_path.is_dir():
        continue
    matrix = matrix_path.name
    plot_dir = Path("plots") / matrix

    plot_time_line(matrix, df_time_summary, plot_dir, schedule="static", chunksize=100)
    plot_cache_misses_bar(matrix, df_perf_summary, plot_dir)
    plot_scheduling_comparison(matrix, df_time_summary, plot_dir)
    print(f"Generated 3 plots for {matrix}")