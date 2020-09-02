#!/bin/bash

script_directory=$(dirname ${BASH_SOURCE[0]})

config=$1
outdir=$2
workdir=${3:-$USERWORK/ellipsis_work}

mkdir -p ${outdir}/config_files
cp ${script_directory}/main.nf ${outdir}/config_files
cp ${config} ${outdir}/config_files

module load Java/11.0.2
nextflow run ${script_directory}/main.nf -c ${config} --out_dir=${outdir} -work-dir ${workdir} -resume
module unload Java/11.0.2
