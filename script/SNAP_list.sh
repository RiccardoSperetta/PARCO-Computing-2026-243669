#!/bin/bash

# Shared graphs list, coming from SNAP - value:
#   1 = assumed undirected, for each a->b edge, also b->a exists
#   0 = assumed directed OR undirected but without the b->a edge
declare -A SNAP_GRAPHS=(
    ["com-dblp.ungraph"]=0
    ["as-skitter"]=0
    ["com-orkut.ungraph"]=0
    ["com-friendster.ungraph"]=0
)