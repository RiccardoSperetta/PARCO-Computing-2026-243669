#!/bin/bash

# Download and extract SNAP graph dataset
# Usage: ./download_snap.sh <dataset_name>
# Example: ./download_snap.sh ca-GrQc
set -e

if [ $# -eq 0 ]; then
    echo "Usage: $0 <dataset_name>"
    echo "Example: $0 ca-GrQc"
    exit 1
fi

DATASET=$1
URL="https://snap.stanford.edu/data/${DATASET}.txt.gz"

echo "Downloading ${DATASET}..."
wget $URL
# Create the target directory if it doesn't exist
mkdir -p data/raw

if [ $? -ne 0 ]; then
    echo "Error: Download failed"
    exit 1
fi

echo "Extracting ${DATASET}.txt.gz..."
gunzip ${DATASET}.txt.gz

# Move the downloaded file to the target directory
mv ${DATASET}.txt data/raw/


if [ $? -eq 0 ]; then
    echo "Success! File extracted to ${DATASET}.txt"
else
    echo "Error: Extraction failed"
    exit 1
fi
