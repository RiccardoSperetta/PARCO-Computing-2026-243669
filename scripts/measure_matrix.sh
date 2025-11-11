#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <matrix_name>"
    exit 1
fi

MATRIX=$1
PROCESSED_DIR="data/processed/${MATRIX}"
RESULTS_BASE="results/${MATRIX}"

# Check if processed data exists
if [ ! -f "${PROCESSED_DIR}/rowptr.bin" ] || \
   [ ! -f "${PROCESSED_DIR}/col.bin" ] || \
   [ ! -f "${PROCESSED_DIR}/val.bin" ]; then
    echo "[ERROR] Processed data not found for $MATRIX"
    exit 1
fi

echo "=== Benchmarking $MATRIX ==="


# BENCHMARK SERIAL CODE: ===============================================
OPT_LEVELS=("O0" "O1" "O2" "O3" "Ofast")

for opt in "${OPT_LEVELS[@]}"; do
    echo "- running $MATRIX - sequential with $opt -"
    RESULTS_DIR="${RESULTS_BASE}/${opt}"
    mkdir -p "$RESULTS_DIR"
    
    # Compile
    gcc -Wall -$opt src/SpMV.c -o spmv_${opt}
    
    # === 1. Time measurements (10 runs) ===
    > "${RESULTS_DIR}/time.txt"  # Clear file first
    
    for i in {1..10}; do
        ./spmv_${opt} "$MATRIX" >> "${RESULTS_DIR}/time.txt"  # APPEND
    done
    
    # === 2. Perf cache misses (10 runs) ===
    > "${RESULTS_DIR}/perf.txt"  # Clear file first
    
    for i in {1..10}; do
        perf stat -e cache-misses,cache-references \
            ./spmv_${opt} "$MATRIX" 2>&1 | \
            grep -E "cache-misses|cache-references" >> "${RESULTS_DIR}/perf.txt"
    done
    
done


# BENCHMARK PARALLEL CODE: ===============================================
# binding each thread to a core, aiming to reduce cache misses
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

# for: number of threads    (2, ..., double the amount of cores)
# = testing for hyperthreading effects
NUM_CORES=$(nproc)
SCHEDULES=("static", "dynamic", "guided")
threads=2

while [ $threads -le $((2 * NUM_CORES)) ]; do
    export OMP_NUM_THREADS=$threads

    # no schedule specified:
    echo "- running $MATRIX - parallel with $threads - default schedule -"
    RESULTS_DIR="${RESULTS_BASE}/th${threads}/default"
    mkdir -p "$RESULTS_DIR"
    
    # Compile
    gcc -Wall -O3 src/SpMV.c -o spmv_${threads}_default -fopenmp
    
    # === 1. Time measurements (10 runs) ===
    > "${RESULTS_DIR}/time.txt"  # Clear file first
    
    for i in {1..10}; do
        ./spmv_${threads}_default "$MATRIX" >> "${RESULTS_DIR}/time.txt"  # APPEND
    done
    
    # === 2. Perf cache misses (10 runs) ===
    > "${RESULTS_DIR}/perf.txt"  # Clear file first
    
    for i in {1..10}; do
        perf stat -e cache-misses,cache-references \
            ./spmv_${threads}_default "$MATRIX" 2>&1 | \
            grep -E "cache-misses|cache-references" >> "${RESULTS_DIR}/perf.txt"
    done


    # changing schedules:



    threads=$((threads * 2))
done



# Cleanup executables
rm -f spmv_*

echo "=== Benchmarking complete for $MATRIX ==="
echo "Results saved in: results/$MATRIX/"

