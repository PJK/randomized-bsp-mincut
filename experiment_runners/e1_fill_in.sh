#!/usr/bin/env bash

reps=$1

cores_lo=( 2 4 8 16 )

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh


out96k=~/evaluation_outputs/e1/96k_extra

mkdir -p ${out96k}

# Fucking 2011 bash cannot pass arrays....
for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	for cores in "${cores_lo[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e1/96k 1:00:00 ${cores} ${out96k} ${seed}
	done
done