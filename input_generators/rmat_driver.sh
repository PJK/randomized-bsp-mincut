#!/usr/bin/env bash

# Usage: bin n m

echo "RMAT driver, `date`"
echo $2 $3
$1 -nVertices $2 -nEdges $3 -noEdgeToSelf -noDuplicateEdges >> /dev/null
gsed -r 's/^([0-9]+)[ \t]*([0-9]+)/\1 \2 1/' out.txt