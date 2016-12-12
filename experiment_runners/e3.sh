#!/usr/bin/env bash

reps=$1

cores_mid=( 72 144 216 252 288 360 396 432 504 576 648 )

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh

out256k=~/evaluation_outputs/e3/256k
out8M=~/evaluation_outputs/e3/8M
out16M=~/evaluation_outputs/e3/16M

mkdir -p ${out256k}
mkdir -p ${out8M}
mkdir -p ${out16M}

for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	for cores in "${cores_mid[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e3/256k 2:00 ${cores} ${out256k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e3/8M 2:00 ${cores} ${out8M} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e3/16M 4:00 ${cores} ${out16M} ${seed}
	done
done
