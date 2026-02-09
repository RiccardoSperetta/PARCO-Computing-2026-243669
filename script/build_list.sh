#!/bin/bash
set -e

source "$(dirname "$0")/SNAP_list.sh"

echo "[ Building CSR for all graphs ]"

for file in data/raw/*; do
    graph=$(basename "$file" .txt)
    csr_file="data/csr/${graph}.bin"
    directed=0 #assuming default edge lists need to be doubled (for each a->b add also b->a)
    shuffle=1 #can be choosen, usually doesn't impact on performance significantly
    csr_file_shuffled="data/csr/${graph}_shuffled.bin"

    if [ ! -f "${file}" ]; then
        echo "Error: ${file} not found"
        continue
    fi
    if [ "${shuffle}" -eq 1 ] && [ -f "${csr_file_shuffled}" ]; then
        echo "Shuffled CSR for ${file} has already been built"
        continue
    elif [ "${shuffle}" -eq 0 ] && [ -f "${csr_file}" ]; then
        echo "CSR for ${file} has already been built"
        continue
    fi

    if [[ -v SNAP_GRAPHS["${graph}"] ]]; then
        directed=${SNAP_GRAPHS["${graph}"]}
    fi

    ./csrMaker.out ${file} ${directed} ${shuffle}

    if [ $? -eq 0 ]; then
        if [ "${shuffle}" -eq 1 ]; then
            echo "Success: CSR built in ${csr_file_shuffled}"
        else
            echo "Success: CSR built in ${csr_file}"
        fi
    else
        echo "Error: CSR construction of ${file} failed"
        exit 1
    fi

    # setup for results folder
    if [[ "${graph}" == kronecker* ]]; then
        mkdir -p "results/weak_scaling"
    else
        if [ "${shuffle}" -eq 1 ]; then
            mkdir -p "results/${graph}_shuffled"
        else
            mkdir -p "results/${graph}"
        fi
    fi
    

done

echo "[ DONE ]"