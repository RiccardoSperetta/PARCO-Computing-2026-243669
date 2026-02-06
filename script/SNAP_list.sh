#!/bin/bash

# Shared graphs list, coming from SNAP - value:
#   1 = assumed undirected, for each a->b edge, also b->a exists
#   0 = assumed directed OR undirected but without the b->a edge => will be turned into undirected graph
declare -A SNAP_GRAPHS=(
    ["com-dblp"]=0
    ["as-skitter"]=0
    ["com-orkut"]=0
    ["ca-AstroPh"]=1
)