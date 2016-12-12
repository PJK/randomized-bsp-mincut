#!/usr/bin/env bash

function run() {
	local input=$1
	local seed=$2
	local cores=$3
	local time=$4
	local output=$5

	echo "#!/usr/bin/env bash
srun --cpu_bind=rank $SCRATCH/square_root 0.95 $input $seed
	 " | sbatch -n ${cores} \
		--time=${time} \
		--ntasks-per-node=36 \
		--output=${output}_%j.out
}

# TODO node tiling for lwo concurrency

function dir_contents() {
	echo $(ls $1/*.in)
}

function process_dir() {
	local dir=$1
	local time=$2
	local cores=$3
	local output_dir=$4
	local seed=$5

	for input in $(dir_contents ${dir}); do
		run ${input} ${seed} ${cores} ${time} ${output_dir}/$(basename ${input})_${cores}
	done
}