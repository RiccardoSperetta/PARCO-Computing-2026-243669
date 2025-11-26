import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# local python scripts imports
from loaders import load_all_data
from summaries import add_summary_stats

PROJECT_ROOT = Path(".")                     
RESULTS_DIR = PROJECT_ROOT / "results"
PLOTS_DIR = PROJECT_ROOT / "plots"


### --------------------------------------------------------------------
# PLOTTING FUNCTIONS
### --------------------------------------------------------------------

### --------------------------------------------------------------------
# shows the 90% percentile of time spent by sequential code(all optimizations) 
# + parallel code with different number of threads(powers of 2)
# chunksize is fixed at 100, but all 3 schedules are being shown
### --------------------------------------------------------------------
def plot_time_line(matrix: str, df_time_summary: pd.DataFrame, plot_dir: Path):
    plot_dir.mkdir(parents=True, exist_ok=True)

    # === Sequential ===
    seq = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["threads"] == 1)
    ].copy()
    seq["xpos"] = seq["opt_level"].map({"O0":0, "O1":1, "O2":2, "O3":3, "Ofast":4})

    # === Parallel (chunksize=100 only) ===
    par = df_time_summary[
        (df_time_summary["matrix"] == matrix) &
        (df_time_summary["chunksize"] == 100) &
        (df_time_summary["threads"] > 1)
    ].copy()
    thread_to_x = {2:5, 4:6, 8:7, 16:8, 32:9, 64:10, 128:11}
    par["xpos"] = par["threads"].map(thread_to_x)

    # === Plot ===
    plt.figure(figsize=(11.5, 5.5))

    # 1. Sequential line (solid, strong blue)
    plt.plot(seq["xpos"], seq["time_ms_p90"],
             'o-', color="#1f77b4", linewidth=3, markersize=8, markerfacecolor="white",
             label="sequential", zorder=10)

    # 2. Parallel lines — slightly transparent + connect from Ofast
    colors = {"static": "#d62728", "dynamic": "#2ca02c", "guided": "#ff7f0e"}
    for sched in ["static", "dynamic", "guided"]:
        sub = par[par["schedule"] == sched]
        if sub.empty:
            continue

        # Get Ofast value to connect from
        ofast_y = seq[seq["opt_level"] == "Ofast"]["time_ms_p90"].iloc[0]

        # Full x and y including the connection from Ofast
        x_full = [4] + sub["xpos"].tolist()        # start from Ofast x=4
        y_full = [ofast_y] + sub["time_ms_p90"].tolist()

        plt.plot(x_full, y_full,
                 'o-', color=colors[sched], linewidth=2.5, markersize=7,
                 alpha=0.85, markerfacecolor="white",
                 label=f"{sched}")

    # === Axis formatting ===
    labels = ["O0","O1","O2","O3","Ofast", "2","4","8","16","32","64","128"]
    plt.xticks(range(12), labels, rotation=45, ha="right")

    plt.ylabel("Execution time (ms) — 90th percentile")
    plt.title(f"Execution time scaling — {matrix}\n(chunksize=100 for parallel runs)")
    plt.grid(True, alpha=0.3)
    plt.xlim(-0.5, 11.5)
    plt.ylim(0, None)
    plt.legend(frameon=False)

    plt.savefig(plot_dir / "time_line_all_schedules.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_speedup_all_matrices(df_time_summary: pd.DataFrame, plot_dir: Path):
    plot_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10.5, 6.5))

    # Get all matrices that have data
    matrices = df_time_summary["matrix"].unique()
    colors = plt.cm.tab10.colors 
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    for i, matrix in enumerate(sorted(matrices)):
        # Baseline = best sequential for this matrix
        baseline = df_time_summary[
            (df_time_summary["matrix"] == matrix) &
            (df_time_summary["threads"] == 1)
        ]["time_ms_p90"].min()

        data = df_time_summary[
            (df_time_summary["matrix"] == matrix) &
            (df_time_summary["schedule"] == "static") &
            (df_time_summary["chunksize"] == 100)
        ].copy()

        if data.empty:
            continue

        speedup = baseline / data["time_ms_p90"]
        threads = data["threads"]

        plt.plot(threads, speedup, 
                 marker=markers[i % len(markers)], 
                 color=colors[i % 10], 
                 linewidth=2.8, markersize=8,
                 label=matrix.replace("_", " "))

    # Ideal line
    threads = [1, 2, 4, 8, 16, 32, 64, 128]
    plt.plot(threads, threads, '--', color='black', linewidth=2, label="Ideal linear")

    plt.xlabel("Number of threads", fontsize=13)
    plt.ylabel("Speedup (vs best sequential)", fontsize=13)
    plt.title("Strong scaling — all matrices\n(OpenMP static, chunksize=100)", fontsize=15)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)
    plt.grid(True, alpha=0.3)
    plt.xscale("log", base=2)
    plt.xticks(threads, threads)
    plt.xlim(1, 128)
    plt.ylim(0, 70)

    plt.tight_layout()
    plt.savefig(plot_dir / "speedup_all_matrices.png", dpi=300, bbox_inches="tight")
    plt.close()

### --------------------------------------------------------------------
# Keeps now only O3 sequential and the parallel static_100 runs
# 2. For every thread count makes two bars side-by-side: orange = L1 miss %, green = LLC miss %
# 3. Height of each bar = 90th percentile (again, realistic worst case)
# 4. Saves one clean bar chart
### --------------------------------------------------------------------
def plot_cache_misses_bar(matrix: str, df_perf_summary: pd.DataFrame, plot_dir: Path):
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
            width, label="L1 miss %", color="#fbff0e", edgecolor="black")
    plt.bar([i + width/2 for i in x], data["LLC_miss_percent_p90"],
            width, label="LLC miss %", color="#a02c77", edgecolor="black")

    plt.xticks(x, data["config"], rotation=45)
    plt.ylabel("Cache miss rate — 90th percentile (%)")
    plt.title(f"Cache miss rates — {matrix}\n(static schedule, chunksize=100)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(plot_dir / "cache_misses.png", dpi=300)
    plt.close()


### --------------------------------------------------------------------
# shows the 90% percentile of time spent by sequential code(all optimizations) 
# + parallel code with different number of threads(powers of 2)
# chunksize is fixed at 100, but all 3 schedules are being shown
### --------------------------------------------------------------------
def plot_scheduling_comparison(matrix: str, df_time_summary: pd.DataFrame, plot_dir: Path):
    threads_list = [8, 16, 32, 64, 128]

    plot_dir.mkdir(parents=True, exist_ok=True)

    # Kept only parallel runs with ≥8 threads for simplicity
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
    colors = ["#d62728", "#2ca02c", "#ff7f0e"]

    for ax, chunksize in zip(axes, chunksizes):
        subset = data[data["chunksize"] == chunksize]

        for i, sched in enumerate(schedules):
            sched_data = subset[subset["schedule"] == sched]
            x_pos = [threads.index(t) + i*bar_width for t in sched_data["threads"]]
            ax.bar(x_pos,
                   sched_data["time_ms_p90"],
                   width=bar_width,
                   label=sched if ax == axes[0] else "",
                   color=colors[i],
                   edgecolor="black")

        ax.set_title(f"chunksize = {chunksize}")
        ax.set_xlabel("Threads")
        ax.set_xticks([threads.index(t) + bar_width for t in threads])
        ax.set_xticklabels(threads)

    axes[0].set_ylabel("Time (ms) — 90th percentile")
    fig.suptitle(f"Scheduling comparison — {matrix}", fontsize=16)
    fig.legend(schedules, loc="upper left", ncol=3)
    plt.tight_layout()
    plt.savefig(plot_dir / "scheduling_comparison.png", dpi=300)
    plt.close()


### --------------------------------------------------------------------
# "MAIN"
### --------------------------------------------------------------------

# ----------------------------------------------------------------------
# From the loaders.py "load_all_data" function I retrieve data frames loaded 
# with all the info stored in the results folder
# ----------------------------------------------------------------------
[df_time, df_perf] = load_all_data(RESULTS_DIR)


# ====================== SUMMARIES WITH 90TH PERCENTILE ======================
# now the previous data frames are summarised, showing their mean and 90% percentile
df_time_summary = add_summary_stats(df_time, ["time_ms"])
df_perf_summary = add_summary_stats(df_perf, ["L1_miss_percent", "LLC_miss_percent"])

PLOTS_DIR.mkdir(parents=True, exist_ok=True)
df_time_summary.to_csv(PLOTS_DIR / "time_summary.csv", index=False)
df_perf_summary.to_csv(PLOTS_DIR / "perf_summary.csv", index=False)

# general speed up obtained and confronted with all matrices on static schedule and chunk-size = 100
plot_speedup_all_matrices(df_time_summary, Path("plots"))

# iterative call on all matrices to build plots:
for matrix_path in RESULTS_DIR.iterdir():
    if not matrix_path.is_dir():
        continue
    matrix = matrix_path.name
    plot_dir = Path("plots") / matrix

    plot_time_line(matrix, df_time_summary, plot_dir)
    plot_cache_misses_bar(matrix, df_perf_summary, plot_dir)
    plot_scheduling_comparison(matrix, df_time_summary, plot_dir)
    print(f"Generated 3 plots for {matrix}")