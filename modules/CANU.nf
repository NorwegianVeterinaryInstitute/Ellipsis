process CANU_NANOPORE {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Canu"

        publishDir "${params.out_dir}/", pattern: "", mode: "copy"
        publishDir "${params.out_dir}/", pattern: "", mode: "copy"

        tag "$datasetID"

        input:
        tuple val(datasetID), file(longreads)

        output:
        file("*")

        """
	canu -correct -p $datasetID genomeSize=$params.genomesize -useGrid=false -nanopore-raw $longreads &> ${datasetID}_canu.log
        """
}

process CANU_PACBIO {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Canu"

        publishDir "${params.out_dir}/", pattern: "", mode: "copy"
        publishDir "${params.out_dir}/", pattern: "", mode: "copy"

        tag "$datasetID"

        input:
        tuple val(datasetID), file(longreads)

        output:
        file("*")

        """
        canu -correct -p $datasetID genomeSize=$params.genomesize -useGrid=false -pacbio-raw $longreads &> ${datasetID}_canu.log
        """
}

