process REPORT {
        module 'R/4.0.0-foss-2020a'

        publishDir "${params.out_dir}/reports", pattern: "total_report.txt", mode: "copy"
        publishDir "${params.out_dir}/reports", pattern: "summary_report.txt", mode: "copy"
        publishDir "${params.out_dir}/reports", pattern: "resfinder_report.txt", mode: "copy"
        publishDir "${params.out_dir}/reports", pattern: "virfinder_report.txt", mode: "copy"
        publishDir "${params.out_dir}/reports", pattern: "plasmidfinder_report.txt", mode: "copy"

        input:
        file("*")
        val(run_ariba)

        output:
        file("*")

        """
        Rscript $baseDir/bin/collate_data.R $run_ariba
        """
}

process REPORT_PLASMAP {
        module 'R/4.0.0-foss-2020a'

        publishDir "${params.out_dir}/plasmap", pattern: "plasmid_mapping_report.txt", mode: "copy"

        input:
        file("*")

        output:
        file("*")

        """
        Rscript $baseDir/bin/collate_plasmap.R
        """
}

