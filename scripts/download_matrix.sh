#!/bin/bash
set -e

# Usage: ./download_matrix.sh <matrix_name> <group>
if [ $# -ne 2 ]; then
    echo "Usage: $0 <matrix_name> <group>"
    echo "Example: $0 1138_bus HB"
    exit 1
fi

BASE_URL="https://suitesparse-collection-website.herokuapp.com/MM"  

MATRIX=$1
GROUP=$2

MTX_FILE="data/raw/${MATRIX}.mtx"

mkdir -p data/raw

# checks if matrix already downloaded:
if [ -f "$MTX_FILE" ]; then
    echo "[SKIP] $MATRIX already exists"
    exit 0
fi

echo "[DOWNLOAD] $MATRIX (group: $GROUP)"

# Download from https://sparse.tamu.edu/:
TEMP_TAR="data/raw/temp_${MATRIX}.tar.gz"
if ! wget -q -O "$TEMP_TAR" "${BASE_URL}/${GROUP}/${MATRIX}.tar.gz"; then
    echo "  [ERROR] Download failed"
    rm -f "$TEMP_TAR"
    exit 1
fi

# Extract to temp directory
TEMP_DIR="data/raw/temp_${MATRIX}"
mkdir -p "$TEMP_DIR"
tar -xzf "$TEMP_TAR" -C "$TEMP_DIR"

# Find .mtx file
MTX_FOUND=$(find "$TEMP_DIR" -name "*.mtx" -o -name "*.MTX" | head -1)

# If found, places it into data/raw/${MATRIX} = $MTX_FILE
if [ -n "$MTX_FOUND" ]; then
    mv "$MTX_FOUND" "$MTX_FILE"
    echo "  [OK]"
else
    echo "  [ERROR] No .mtx file found"
    echo "  Contents:"
    ls -la "$TEMP_DIR"
fi

# Cleanup
rm -rf "$TEMP_TAR" "$TEMP_DIR"