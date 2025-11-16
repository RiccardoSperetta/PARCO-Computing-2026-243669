#!/bin/bash
set -e

source "$(dirname "$0")/matrix_list.sh"

echo "[ Processing matrices ]"

if [ ! -f "src/matrix_processing.c" ]; then
    echo "[ERROR] process_matrix program not found"
    exit 1
fi

for matrix in "${!MATRIX_INFO[@]}"; do
    ./scripts/pre_process_matrix.sh "$matrix"
done

echo "[ DONE ]"
