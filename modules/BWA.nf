process BWA_MEM_SAM {
	module 'BWA/0.7.17-foss-2018b'
	module 'SAMtools/1.9-foss-2018b'

        publishDir "${params.out_dir}/plasmap/", pattern: "*.bam", mode: "copy"

        input:
        tuple val(datasetID), file(R1), file(R2), val(ref_id), file(ref)

        output:
        file("*")
	tuple val(datasetID), path("*.bam"), val(ref_id), emit: bam_ch

        script:
        """
	bwa index $ref
	bwa mem -t $task.cpus $ref $R1 $R2 | samtools sort -o ${ref_id}_${datasetID}_mapped_sorted.bam
        """
}
