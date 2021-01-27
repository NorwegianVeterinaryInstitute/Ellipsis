process CANU {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Canu"

	publishDir "${params.out_dir}/results/canu", pattern: "*report", mode: "copy"
	publishDir "${params.out_dir}/results/canu", pattern: "*log", mode: "copy"

        tag "$datasetID"
	label 'hybrid'

        input:
	val type
        tuple val(datasetID), file(longreads)

        output:
        file("*")
	tuple val(datasetID), path {"*correctedReads.fasta.gz"}, emit: canu_output

	script:
	if (type == "nanopore") {
        	"""
		canu -correct -p $datasetID genomeSize=$params.genomesize -useGrid=false -nanopore-raw $longreads &> ${datasetID}_canu.log
        	"""
	}
	if (type == "pacbio") {
		"""
		canu -correct -p $datasetID genomeSize=$params.genomesize -useGrid=false -pacbio-raw $longreads &> ${datasetID}_canu.log
		"""
	}
}
