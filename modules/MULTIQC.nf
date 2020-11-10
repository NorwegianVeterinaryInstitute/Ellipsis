process MULTIQC {
        module "FastQC/0.11.9-Java-11 "

        publishDir "${params.out_dir}/fastqc/", pattern: "*", mode: "copy"

        tag {$datasetID}

        input:
        tuple val(datasetID), file(reads)

        output:
        file("*")
        path "$datasetID", emit: fastqc_reports


        """
        mkdir $datasetID
        fastqc $reads -o $datasetID -t $task.cpus
        """
}

