#!/usr/bin/env bash

reps=$1

cores_lo=( 1 4 8 16 ) # These need a wrapper.....
cores_mid=( 36 72 108 144 288 432 )
cores_hi=( 576 720 864 1008 1152 1296 )

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh

mkdir -p ~/evaluation_outputs/e2

for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	for cores in "${cores_mid[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e2 1:00 ${cores} ~/evaluation_outputs/e2 ${seed}
	done

	for cores in "${cores_hi[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e2 1:00 ${cores} ~/evaluation_outputs/e2 ${seed}
	done
done