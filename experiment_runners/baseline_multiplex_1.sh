#!/usr/bin/env bash

seed=$RANDOM


function run_two() {
	local input1=$1
	local input2=$2
	local time=$3
	local output=$4


	echo "#!/usr/bin/env bash
$(for i in $(seq 1 18); do
	echo "$SCRATCH/karger_stein $input1 $(($seed + $i)) &"
done)

$(for i in $(seq 1 18); do
	echo "$SCRATCH/karger_stein $input2 $(($seed + 0)) &"
done)

wait
		 " | sbatch --time=${time} \
			-N 1 \
			--output=${output}/basline_%j.out
}

# hi memory
function run_one() {
	local input1=$1
	local time=$2
	local output=$3


	echo "#!/usr/bin/env bash
$(for i in $(seq 1 6); do
	echo "$SCRATCH/karger_stein $input1 $(($seed + $i)) &"
done)
wait
		 " | sbatch --time=${time} \
			-N 1 \
			--output=${output}/basline_%j.out
}

baselines1=~/evaluation_outputs/e1/baselines
mkdir -p ${baselines1}
#run_two $SCRATCH/evaluation_inputs/e1/32k/ws_32k_16.in $SCRATCH/evaluation_inputs/e1/48k/ws_48k_16.in 3:00 ${baselines1}
#run_two $SCRATCH/evaluation_inputs/e1/64k/ws_64k_16.in $SCRATCH/evaluation_inputs/e1/96k/ws_96k_16.in 5:00 ${baselines1}


baselines3=~/evaluation_outputs/e3/baselines
mkdir -p ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/16M/ba_16k_16M.in 10:00 ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/16M/ws_16k_16M.in 10:00 ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/256k/ba_16k_256k.in 10:00 ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/256k/ws_16k_265k.in 10:00 ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/8M/ba_16k_8M.in 10:00 ${baselines3}
run_one $SCRATCH/evaluation_inputs/e3/8M/ws_16k_8M.in 10:00 ${baselines3}

baselines4=~/evaluation_outputs/e4/baselines
mkdir -p ${baselines4}
#run_one $SCRATCH/evaluation_inputs/e4/rmat_16k_1024.in 2:00 ${baselines4}
#run_one $SCRATCH/evaluation_inputs/e4/rmat_16k_2048.in 3:00 ${baselines4}
#run_one $SCRATCH/evaluation_inputs/e4/rmat_16k_4096.in 4:00 ${baselines4}

