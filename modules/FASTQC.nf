process FASTQC {
	module="FastQC/0.11.8-Java-1.8"

        publishDir "${params.out_dir}/fastqc/", pattern: "$datasetID", mode: "copy"
	
	tag "$datasetID"

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
	path "$datasetID", emit: fastqc_reports


        """
	mkdir $datasetID
	fastqc $R1 $R2 -o $datasetID -t $task.cpus
        """
}

