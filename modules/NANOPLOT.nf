process NANOPLOT {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/nanoplot"

        publishDir "${params.out_dir}/results/nanoplot/", pattern: "$datasetID", mode: "copy"
	
	tag "$datasetID"

        input:
        tuple val(datasetID), file(reads)

        output:
        file("*")

        """
	mkdir $datasetID
	NanoPlot -t $task.cpus --fastq $reads --N50 -o $datasetID
        """
}

