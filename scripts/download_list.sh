#!/bin/bash
set -e

# Load the shared matrix list
source "$(dirname "$0")/matrix_list.sh"

echo "[ Downloading matrices ]"

for matrix in "${!MATRIX_INFO[@]}"; do
    GROUP="${MATRIX_INFO[$matrix]}"
    
    # Call the single-matrix download script
    ./scripts/download_matrix.sh "$matrix" "$GROUP"
done

echo "[ DONE ]"