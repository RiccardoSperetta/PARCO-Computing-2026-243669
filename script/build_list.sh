#!/bin/bash

set -e

echo "[ Building CSR for all graphs ]"

for file in data/raw/*; do
    graph=$(basename "$file" .txt)
    csr_file="data/csr/${graph}.bin"

    if [ ! -f "${file}" ]; then
        echo "Error: ${file} not found"
        exit 1
    elif [ -f "${csr_file}" ]; then
        echo "CSR for ${file} has already been built"
        continue
    fi

    ./a.out ${file} 1 1 # usually assuming undirected graphs - not shuffling

    if [ $? -eq 0 ]; then
        echo "Success: CSR built in ${csr_file}"
    else
        echo "Error: CSR construction of ${file} failed"
        exit 1
    fi
done

echo "[ DONE ]"
