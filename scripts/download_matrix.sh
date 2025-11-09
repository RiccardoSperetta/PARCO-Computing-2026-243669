#!/bin/bash
set -e

declare -A MATRIX_INFO=(
    ["1138_bus"]="HB"
    ["crankseg_2"]="GHS_psdef"
    # up to five!
)

BASE_URL="https://suitesparse-collection-website.herokuapp.com/MM"
mkdir -p data/raw

echo "=== Downloading matrices ==="

for matrix in "${!MATRIX_INFO[@]}"; do
    GROUP="${MATRIX_INFO[$matrix]}"
    MTX_FILE="data/raw/${matrix}.mtx"
    
    if [ -f "$MTX_FILE" ]; then
        echo "[SKIP] $matrix"
        continue
    fi
    
    echo "[DOWNLOAD] $matrix"
    
    TEMP_TAR="data/raw/temp_${matrix}.tar.gz"
    if ! wget -q -O "$TEMP_TAR" "${BASE_URL}/${GROUP}/${matrix}.tar.gz"; then
        echo "  [ERROR] Download failed"
        rm -f "$TEMP_TAR"
        continue
    fi
    
    # Extract to temp directory
    TEMP_DIR="data/raw/temp_${matrix}"
    mkdir -p "$TEMP_DIR"
    tar -xzf "$TEMP_TAR" -C "$TEMP_DIR"
    
    # Find .mtx file
    MTX_FOUND=$(find "$TEMP_DIR" -name "*.mtx" -o -name "*.MTX" | head -1)
    
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
done

echo "=== Done ==="