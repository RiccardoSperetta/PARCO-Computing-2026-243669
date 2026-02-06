#!/bin/bash
set -e

for dir in results/*; do
    rm -rf "$dir"/*

done

echo "[ DONE ]"