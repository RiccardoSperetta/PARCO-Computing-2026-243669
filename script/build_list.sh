#!/bin/bash
set -e

source "$(dirname "$0")/SNAP_list.sh"

echo "[ Building CSR for all graphs ]"

for file in data/raw/*; do
    graph=$(basename "$file" .txt)
    csr_file="data/csr/${graph}.bin"
    undirected=0

    if [ ! -f "${file}" ]; then
        echo "Error: ${file} not found"
        continue
    elif [ -f "${csr_file}" ]; then
        echo "CSR for ${file} has already been built"
        continue
    fi

    if [[ -v SNAP_GRAPHS["${graph}"] ]]; then
        undirected=${SNAP_GRAPHS["${graph}"]}
    fi

    ./csrMaker.out ${file} ${undirected} 0 # - not shuffling

    if [ $? -eq 0 ]; then
        echo "Success: CSR built in ${csr_file}"
    else
        echo "Error: CSR construction of ${file} failed"
        exit 1
    fi

    # setup for results folder
    if [[ "${graph}" == kronecker* ]]; then
        mkdir -p "results/weak_scaling"
    else
        mkdir -p "results/${graph}"
    fi
    

done

echo "[ DONE ]"