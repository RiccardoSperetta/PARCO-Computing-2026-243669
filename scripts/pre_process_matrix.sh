#!/bin/bash
set -e

# Usage: ./process_matrix.sh <matrix_name>

if [ $# -ne 1 ]; then
    echo "Usage: $0 <matrix_name>"
    echo "Example: $0 1138_bus"
    exit 1
fi

MATRIX=$1
RAW_FILE="data/raw/${MATRIX}.mtx"
PROCESSED_DIR="data/processed/${MATRIX}"

if [ ! -f "$RAW_FILE" ]; then
    echo "[ERROR] Raw file not found: $RAW_FILE, must download it first"
    exit 1
fi

if [ -f "${PROCESSED_DIR}/rowptr.bin" ] && \
   [ -f "${PROCESSED_DIR}/col.bin" ] && \
   [ -f "${PROCESSED_DIR}/val.bin" ]; then
    echo "[SKIP] $MATRIX - already processed"
    exit 0
fi

echo "[PROCESS] $MATRIX"

mkdir -p "$PROCESSED_DIR"

if ./src/matrix_processing.out "${MATRIX}.mtx"; then
    echo "  [OK]"
else
    echo "  [ERROR] Processing failed"
    rm -f "${PROCESSED_DIR}"/*.bin
    exit 1
fi