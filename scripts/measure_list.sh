#!/bin/bash
set -e

source "$(dirname "$0")/matrix_list.sh"

echo "[ Measuring matrices ]"

#pointless to run if no matrix has been processed yet
if [ ! -d "data/processed" ]; then
    echo "[ERROR] processed matrix folder doesn't exist yet"
    exit 1
fi

#pointless to run if there's no code to run
if [ ! -f "src/SpMV.c" ]; then
    echo "[ERROR] missing C code to compile and run"
    exit 1
fi

for matrix in "${!MATRIX_INFO[@]}"; do
    ./scripts/measure_matrix.sh "$matrix"
done

echo "[ Measuring complete ]"