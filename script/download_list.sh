#!/bin/bash

set -e
# Load the shared graphs list
source "$(dirname "$0")/SNAP_list.sh"

echo "[ Downloading graphs ]"

for graph in "${SNAP_GRAPHS[@]}"; do
    if [ -f "data/raw/${graph}.txt" ]; then
        echo "File data/raw/${graph}.txt already exists."
        continue
    fi

    URL="https://snap.stanford.edu/data/${graph}.txt.gz"

    echo "Downloading ${graph}..."
    wget $URL
    # Create the target directory if it doesn't exist
    mkdir -p data/raw

    if [ $? -ne 0 ]; then
        echo "Error: Download failed"
        exit 1
    fi

    gunzip ${graph}.txt.gz

    # Move the downloaded file to the target directory
    mv ${graph}.txt data/raw/


    if [ $? -eq 0 ]; then
        echo "Success: File extracted to ${graph}.txt"
    else
        echo "Error: Extraction failed"
        exit 1
    fi

done

echo "[ DONE ]"
