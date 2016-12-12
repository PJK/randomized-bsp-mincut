#!/usr/bin/env bash

reps=$1

cores_lo=( 1 )
cores_mid=( 36 72 108 144 288 432 )
cores_hi=( 576 720 864 1008 1152 1296 1440 )

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh


out32k=~/evaluation_outputs/e1/32k
out48k=~/evaluation_outputs/e1/48k
out64k=~/evaluation_outputs/e1/64k
out96k=~/evaluation_outputs/e1/96k

mkdir -p ${out32k}
mkdir -p ${out48k}
mkdir -p ${out64k}
mkdir -p ${out96k}

for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	for cores in "${cores_lo[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e1/32k 8:00 ${cores} ${out32k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/48k 16:00 ${cores} ${out48k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/64k 40:00 ${cores} ${out64k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/96k 1:00:00 ${cores} ${out96k} ${seed}
	done

	for cores in "${cores_mid[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e1/32k 1:00 ${cores} ${out32k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/48k 2:00 ${cores} ${out48k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/64k 2:00 ${cores} ${out64k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/96k 6:00 ${cores} ${out96k} ${seed}
	done

	for cores in "${cores_hi[@]}"; do
		process_dir $SCRATCH/evaluation_inputs/e1/32k 1:00 ${cores} ${out32k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/48k 1:00 ${cores} ${out48k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/64k 1:00 ${cores} ${out64k} ${seed}
		process_dir $SCRATCH/evaluation_inputs/e1/96k 1:00 ${cores} ${out96k} ${seed}
	done
done