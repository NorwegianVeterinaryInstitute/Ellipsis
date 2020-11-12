process FILTLONG {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Filtlong"

        tag "$datasetID"
        label 'heavy'

        input:
        tuple val(datasetID), file(R1), file(R2), file(longreads)

        output:
        file("*")
        tuple val(datasetID), path {"*filtered.fastq.gz"}, emit: filtered_longreads

        """
	filtlong -1 $R1 -2 $R2 --min_length $params.minlen --keep_percent $params.keep_percent --target_bases $params.target_bases $longreads | gzip > ${datasetID}_filtered.fastq.gz
        """
}
