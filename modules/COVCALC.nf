process COV_CALC {
        module 'BEDTools/2.27.1-foss-2018b'

        publishDir "${params.out_dir}/plasmap/coverage_reports", pattern: "*_genomecov.txt", mode: "copy"

        input:
        tuple val(datasetID), file(bam), val(ref_id)

        output:
        file("*")
	path"*.txt", emit: R_ch

        script:
        """
	bedtools genomecov -ibam $bam -d > ${datasetID}_${ref_id}_genomecov.txt
        """
}
