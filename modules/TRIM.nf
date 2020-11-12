process TRIM {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Trimgalore"

	tag "$datasetID"

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
        tuple val(datasetID), path {"*val_1.fq.gz"}, path {"*val_2.fq.gz"}, emit: trim_reads

        script:
        """
        trim_galore -o . --paired -trim1 $R1 $R2 &> ${datasetID}_trimgalore.log
        """
}

