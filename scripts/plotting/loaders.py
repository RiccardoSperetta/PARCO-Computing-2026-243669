# loaders.py
from pathlib import Path
from typing import List, Dict

import pandas as pd


# Expected values (adjust only if your actual folders differ)
OPT_LEVELS = ["O0", "O1", "O2", "O3", "Ofast"]
THREAD_COUNTS = [2, 4, 8, 16, 32, 64, 128]
SCHEDULES = ["static", "dynamic", "guided"]
CHUNKSIZES = [10, 100, 1000]


# ----------------------------------------------------------------------
# reads time.txt (x values in milliseconds, scientific notation)
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
# reads perf.txt (x lines of "L1_miss% LLC_miss%")
# ----------------------------------------------------------------------
def read_perf(perf_file: Path) -> List[Dict]: 
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


def load_all_data(results_dir: Path):
    # ----------------------------------------------------------------------
    # Storage for all raw measurements
    # ----------------------------------------------------------------------
    time_records: List[Dict] = []
    perf_records: List[Dict] = []

    # ----------------------------------------------------------------------
    # Walk through all matrices: gathering all data in 2 separate data frames, 
    # one for time measurements, one for cache misses
    # ----------------------------------------------------------------------

    for matrix_path in results_dir.iterdir():
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

    # Build final DataFrames from raw data
    return pd.DataFrame(time_records), pd.DataFrame(perf_records) # returned as a tuple
        