import os
import pandas as pd
import matplotlib.pyplot as plt


RESULTS_DIR = "results"
PLOTS_DIR   = "plots"

# -------------------------------------------------------------------
# Utility: ensure a directory exists
# -------------------------------------------------------------------
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


# -------------------------------------------------------------------
# Parse the folder structure to extract metadata
# Example:
# results/matrix1/th8/static_100/time.txt
# → matrix = "matrix1"
# → threads = 8
# → schedule = "static"
# → chunk = 100
# -------------------------------------------------------------------
def parse_metadata(path_parts):
    """
    path_parts example:
    ["matrix1", "th8", "static_100"]

    returns a dict:
    {
        "matrix": "matrix1",
        "threads": 8,
        "schedule": "static",
        "chunk": 100
    }
    """
    matrix = path_parts[0]

    # threads folder is "thX"
    if path_parts[1].startswith("th"):
        threads = int(path_parts[1][2:])
        is_parallel = True
    else:
        # sequential folder
        threads = 1
        is_parallel = False

    # schedule/chunk folder
    sched_folder = path_parts[2]

    if "_" in sched_folder:
        schedule, chunk = sched_folder.split("_")
        chunk = int(chunk)
    else:
        # sequential has only O0/, O3/, etc.
        schedule = sched_folder
        chunk = None

    return {
        "matrix": matrix,
        "threads": threads,
        "schedule": schedule,
        "chunk": chunk,
        "is_parallel": is_parallel,
    }


# -------------------------------------------------------------------
# Read numeric files (time.txt or perf.txt)
# time.txt contains ONE number per line → time in ms or seconds
# perf.txt contains ONE number per line → cache misses or miss rate
# -------------------------------------------------------------------
def read_numeric_file(filepath):
    values = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    values.append(float(line))
                except ValueError:
                    pass  # ignore garbage
    return values


# -------------------------------------------------------------------
# Create a simple plot: metric vs threads
# -------------------------------------------------------------------
def plot_vs_threads(df, metric, matrix_name):
    plt.figure()
    subset = df[df["matrix"] == matrix_name]

    # Group by schedule (static, dynamic, etc.)
    for sched, group in subset.groupby("schedule"):
        plt.plot(group["threads"], group[metric],
                 marker='o', label=sched)

    plt.xlabel("Threads")
    plt.ylabel(metric)
    plt.title(f"{matrix_name}: {metric} vs Threads")
    plt.legend()
    ensure_dir(f"{PLOTS_DIR}/{matrix_name}")
    plt.savefig(f"{PLOTS_DIR}/{matrix_name}/{metric}_vs_threads.png")
    plt.close()


# -------------------------------------------------------------------
# MAIN SCRIPT
# -------------------------------------------------------------------
def main():
    entries = []  # will be a table of all measurements

    for matrix in os.listdir(RESULTS_DIR):
        matrix_path = os.path.join(RESULTS_DIR, matrix)
        if not os.path.isdir(matrix_path):
            continue

        # thX folders or "sequential/"
        for th_folder in os.listdir(matrix_path):
            th_path = os.path.join(matrix_path, th_folder)
            if not os.path.isdir(th_path):
                continue

            # inside → schedule folders like static_100/
            for sched_folder in os.listdir(th_path):
                sched_path = os.path.join(th_path, sched_folder)
                if not os.path.isdir(sched_path):
                    continue

                meta = parse_metadata([matrix, th_folder, sched_folder])

                # Read time.txt
                time_file = os.path.join(sched_path, "time.txt")
                perf_file = os.path.join(sched_path, "perf.txt")

                times = read_numeric_file(time_file) if os.path.exists(time_file) else []
                perf = read_numeric_file(perf_file) if os.path.exists(perf_file) else []

                if len(times) > 0:
                    meta["Time (ms)"] = sum(times) / len(times)
                else:
                    meta["Time (ms)"] = None

                if len(perf) > 0:
                    meta["Cache Misses(%)"] = sum(perf) / len(perf)
                else:
                    meta["Cache Misses(%)"] = None

                entries.append(meta)

    # Convert everything to a DataFrame
    df = pd.DataFrame(entries)

    # Generate plots for each matrix
    for matrix in df["matrix"].unique():
        plot_vs_threads(df, "Time (ms)", matrix)
        plot_vs_threads(df, "Cache Misses(%)", matrix)

    print("Plots generated in:", PLOTS_DIR)


if __name__ == "__main__":
    main()
