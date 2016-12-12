#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function fetch() {
	cat ${DIR}/../evaluation/data/header.csv > ${DIR}/../evaluation/data/$1
	echo >> ${DIR}/../evaluation/data/$1
	ssh kalvodap@ela.cscs.ch "ssh kalvodap@dora.cscs.ch \"cat evaluation_outputs/${2} | sort -V\"" >> ${DIR}/../evaluation/data/$1
}

fetch e1/32k.csv 'e1/32k/*'
fetch e1/48k.csv 'e1/48k/*'
fetch e1/64k.csv 'e1/64k/*'
fetch e1/96k.csv 'e1/96k/*'
fetch e1/96k_extras.csv 'e1/96k_extra/*'
fetch e1/baselines.csv 'e1/baselines/*'

fetch e3/256k.csv 'e3/256k/*'
fetch e3/8M.csv 'e3/8M/*'
fetch e3/16M.csv 'e3/16M/*'
fetch e3/baselines.csv 'e3/baselines/*'

fetch e4/data.csv 'e4/*'
fetch e4/baselines.csv 'e4/baselines/*'

fetch e6/data.csv 'e6/*'