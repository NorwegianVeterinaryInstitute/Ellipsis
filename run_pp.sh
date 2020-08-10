#!/bin/bash
# Run nextflow pipeline
module load Java/11.0.2

cd /cluster/projects/nn9305k/active/hkaspersen/Projects/NEXTFLOW/plasmid_pipeline
nextflow run plasmid_pipeline.nf -c plasmid_pipeline.config -work-dir ${USERWORK}/pp_work -resume

module unload
