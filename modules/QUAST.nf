process QUAST {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/reports/", pattern: "transposed_report.tsv", mode: "copy"

        input:
        file("*")

        output:
        file("*")
        path "transposed_report.tsv", emit: R_quast

        script:
        """
        quast --threads $task.cpus -o . *.fasta
        """
}

