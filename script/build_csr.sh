#!/bin/bash

# Runs csrMaker.c on the specified argument (graph)
set -e

if [ $# -eq 0 ]; then
    echo "Usage: $0 <graph_name>"
    echo "Example: $0 ca-GrQc"
    exit 1
fi

if [ ! -f "data/raw/$1.txt" ]; then
    echo "Error: data/raw/$1.txt not found"
    exit 1
fi

#TODO: make the shuffling a parameter of the script?
./csrMaker data/raw/$1.txt 0 0 # usually assuming undirected graphs 

if [ $? -eq 0 ]; then
    echo "Success: CSR built in data/csr/$1.bin"
else
    echo "Error: CSR construction failed"
    exit 1
fi