#!/bin/bash
set -e

source "$(dirname "$0")/SNAP_list.sh"

echo "[ Building CSR for all graphs ]"

for file in data/raw/*; do
    graph=$(basename "$file" .txt)
    csr_file="data/csr/${graph}.bin"
    directed=0 #assuming default edge lists need to be doubled (for each a->b add also b->a)
    shuffle=0 #can be choosen, shuffling usually produces slightly worse performance

    if [ ! -f "${file}" ]; then
        echo "Error: ${file} not found"
        continue
    fi
    if [ -f "${csr_file}" ]; then
        echo "CSR for ${file} has already been built"
        continue
    fi

    if [[ -v SNAP_GRAPHS["${graph}"] ]]; then
        directed=${SNAP_GRAPHS["${graph}"]}
    fi

    ./csrMaker.out ${file} ${directed} ${shuffle}

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