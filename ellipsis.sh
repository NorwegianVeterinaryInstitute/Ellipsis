#!/bin/bash

script_directory=$(dirname ${BASH_SOURCE[0]})

track=$1
config=$2
outdir=$3
workdir=${4:-$USERWORK/ellipsis_work}

mkdir -p ${outdir}/config_files
cp ${config} ${outdir}/config_files
git --git-dir ${script_directory}/.git branch -v | grep "\*" | awk '{print $2, $3}' > ${outdir}/config_files/ellipsis_version.log

if [[ $track == "main" ]]; then
	cp ${script_directory}/main.nf ${outdir}/config_files

	module load Java/11.0.2
	nextflow run ${script_directory}/main.nf -c ${config} --out_dir=${outdir} -work-dir ${workdir} -resume
	module unload Java/11.0.2; fi

if [[ $track == "plasmap" ]]; then
        cp ${script_directory}/plasmap.nf ${outdir}/config_files

        module load Java/11.0.2
        nextflow run ${script_directory}/plasmap.nf -c ${config} --out_dir=${outdir} -work-dir ${workdir} -resume
        module unload Java/11.0.2; fi
