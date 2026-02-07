#!/bin/bash
set -e

for dir in results/*; do
    rm -rf "$dir"/$1*

done

echo "[ DONE ]"