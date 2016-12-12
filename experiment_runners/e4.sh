#!/usr/bin/env bash

reps=$1
core_counts=( 144 288 432 576 720 864 1008 1152 1296 )

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh


out=~/evaluation_outputs/e4

mkdir -p ${out}

for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	for cores in "${core_counts[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e4 0:55 ${cores} ${out} ${seed}
	done
done