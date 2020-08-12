#!/bin/bash

# input parameters
config=$1

# Run nextflow pipeline
module load Java/11.0.2

nextflow run /cluster/projects/nn9305k/active/hkaspersen/Projects/NEXTFLOW/plasmid_pipeline/plasmid_pipeline.nf -c $config -work-dir ${USERWORK}/pp_work -resume

module unload
