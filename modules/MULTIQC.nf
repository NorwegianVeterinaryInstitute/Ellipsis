process MULTIQC {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/multiqc/", pattern: "*html", mode: "copy"

        input:
        file("*")

        output:
        file("*")

        """
        multiqc *.zip
        """
}

