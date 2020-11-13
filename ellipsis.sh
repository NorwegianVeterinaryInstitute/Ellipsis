#!/bin/bash

script_directory=$(dirname ${BASH_SOURCE[0]})

track=$1
config=$2
outdir=$3
workdir=${4:-$USERWORK/ellipsis_work}

if [[ $track == "main" ]]; then
	mkdir -p ${outdir}/config_files
	cp ${script_directory}/main.nf ${outdir}/config_files
	cp ${config} ${outdir}/config_files

	module load Java/11.0.2
	nextflow run ${script_directory}/main.nf -c ${config} --out_dir=${outdir} -work-dir ${workdir} -resume
	module unload Java/11.0.2; fi

if [[ $track == "plasmap" ]]; then
	mkdir -p ${outdir}/config_files
        cp ${script_directory}/plasmap.nf ${outdir}/config_files
        cp ${config} ${outdir}/config_files

        module load Java/11.0.2
        nextflow run ${script_directory}/plasmap.nf -c ${config} --out_dir=${outdir} -work-dir ${workdir} -resume
        module unload Java/11.0.2; fi
