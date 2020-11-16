process MULTIQC_PRE {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/reports/", pattern: "*html", mode: "copy", saveAs: {"MultiQC_pre_trimming_report.html"}

        input:
        file("*")

        output:
        file("*")

        """
        multiqc *.zip
        """
}

process MULTIQC_POST {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/reports/", pattern: "*html", mode: "copy", saveAs: {"MultiQC_post_trimming_report.html"}

        input:
        file("*")

        output:
        file("*")

        """
        multiqc *.zip
        """
}

