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

# ======================================================================
# BENCHMARK SEQUENTIAL CODE: ===========================================
# ======================================================================
SEQUENTIAL_FLAGS="--std=c11 -Wall src/SpMV.c -o spmv.out"

OPT_LEVELS=("O0" "O1" "O2" "O3" "Ofast")

for opt in "${OPT_LEVELS[@]}"; do
    echo "- running $MATRIX - sequential with $opt -"
    RESULTS_DIR="${RESULTS_BASE}/sequential/${opt}"
    mkdir -p "$RESULTS_DIR"
    
    # Compile
    gcc ${SEQUENTIAL_FLAGS} -$opt 
    
    # === 1. Time measurements (10 runs) ===
    > "${RESULTS_DIR}/time.txt"  # Clear file first
    
    for i in {1..10}; do
        ./spmv.out "$MATRIX" >> "${RESULTS_DIR}/time.txt"  # APPEND
    done
    
    # === 2. Perf cache misses (10 runs) ===
    > "${RESULTS_DIR}/perf.txt"  # Clear file first
    
    for i in {1..10}; do
        perf stat -e cache-misses,cache-references \
            ./spmv.out "$MATRIX" 2>&1 | \
            awk '/cache-misses/ {gsub(/,/, "", $1); misses=$1} 
                /cache-references/ {gsub(/,/, "", $1); refs=$1} 
                END {if (refs > 0) printf "%.4f\n", (misses/refs)*100}' \
            >> "${RESULTS_DIR}/perf.txt"
    done
    
        echo "- done $MATRIX - sequential with $opt -"
done

# ======================================================================
# BENCHMARK PARALLEL CODE: =============================================
# ======================================================================

# binding each thread to a core, aiming to reduce cache misses
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

PARALLEL_FLAGS=${SEQUENTIAL_FLAGS}" -O3 -fopenmp"

# for: number of threads    (2, ..., double the amount of cores)
# = testing for hyperthreading effects
NUM_CORES=$(nproc)
SCHEDULES=("static" "dynamic" "guided")
CHUNKSIZES=("10" "100" "1000")
threads=2 #initial value

while [ $threads -le $((2 * NUM_CORES)) ]; do
    export OMP_NUM_THREADS=$threads

    # for each schedule type
    for sched in ${SCHEDULES[@]}; do
        # and for a selection of chunksizes
        for chunk in ${CHUNKSIZES[@]}; do
            export OMP_SCHEDULE="${sched},${chunk}"

            RESULTS_DIR="${RESULTS_BASE}/th${threads}/${sched}_${chunk}"
            mkdir -p "$RESULTS_DIR"
            
            # Compile
            gcc ${PARALLEL_FLAGS}

            echo "- running $MATRIX - parallel with $threads threads - ${sched} schedule - chunksize = ${chunk} -"            
            # === 1. Time measurements (10 runs) ===
            > "${RESULTS_DIR}/time.txt"  # Clear file first
            
            for i in {1..10}; do
                ./spmv.out "$MATRIX" >> "${RESULTS_DIR}/time.txt"  # APPEND
            done
            
            # === 2. Perf cache misses (10 runs) ===
            > "${RESULTS_DIR}/perf.txt"  # Clear file first
            
            for i in {1..10}; do
                perf stat -e cache-misses,cache-references \
                    ./spmv.out "$MATRIX" 2>&1 | \
                    awk '/cache-misses/ {gsub(/,/, "", $1); misses=$1} 
                        /cache-references/ {gsub(/,/, "", $1); refs=$1} 
                        END {if (refs > 0) printf "%.4f\n", (misses/refs)*100}' \
                    >> "${RESULTS_DIR}/perf.txt"
            done

            echo "- done $MATRIX - parallel with $threads threads - ${sched} schedule - chunksize = ${chunk} -"
            echo ""
        done # chunksizes change done

    done # schedules change done

    threads=$((threads * 2)) #doubling the amount of threads every time, 
done # number of threads change done



# Cleanup executable
rm spmv.out

echo "=== Benchmarking complete for $MATRIX ==="
echo "Results saved in: results/$MATRIX/"

