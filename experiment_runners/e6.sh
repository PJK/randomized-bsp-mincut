#!/usr/bin/env bash

reps=$1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/common.sh


out=~/evaluation_outputs/e6

mkdir -p ${out}

for _rep in $(seq 1 ${reps}); do
	seed=$RANDOM

	run $SCRATCH/evaluation_inputs/e6/ws_32k_16.in ${seed} 144 3:00 ${out}/32
	run $SCRATCH/evaluation_inputs/e6/ws_48k_16.in ${seed} 216 3:00 ${out}/48
	run $SCRATCH/evaluation_inputs/e6/ws_64k_16.in ${seed} 288 3:00 ${out}/64
	run $SCRATCH/evaluation_inputs/e6/ws_80k_16.in ${seed} 360 3:00 ${out}/80
	run $SCRATCH/evaluation_inputs/e6/ws_96k_16.in ${seed} 432 3:00 ${out}/96
done