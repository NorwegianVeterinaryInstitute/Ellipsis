process COV_CALC {
        module 'BEDTools/2.27.1-foss-2018b'

        publishDir "${params.out_dir}/plasmap/", pattern: "*_report.txt", mode: "copy"

        input:
        tuple val(datasetID), file(bam), val(ref_id)

        output:
        file("*")
	path"*mapping_report.txt", emit: R_ch

        script:
        """
	calc_cov_ref.bash $bam $ref_id > ${ref_id}_${datasetID}_mapping_report.txt
        """
}

