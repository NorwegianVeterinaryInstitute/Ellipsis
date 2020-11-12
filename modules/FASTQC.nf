process FASTQC {
	module="FastQC/0.11.9-Java-11"

        publishDir "${params.out_dir}/results/fastqc/", pattern: "$datasetID", mode: "copy"
	
	tag "$datasetID"

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")

        """
	mkdir $datasetID
	fastqc $R1 $R2 -o $datasetID -t $task.cpus
        """
}

