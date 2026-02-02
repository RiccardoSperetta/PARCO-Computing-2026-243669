#!/bin/bash

set -e
# Load the shared graphs list
source "$(dirname "$0")/SNAP_list.sh"

echo "[ Downloading graphs ]"

for graph in "${!SNAP_GRAPHS[@]}"; do
    echo $GRAPH
    
    # Call the download_graph.sh script
    ./script/download_graph.sh "$graph"
done

echo "[ DONE ]"
