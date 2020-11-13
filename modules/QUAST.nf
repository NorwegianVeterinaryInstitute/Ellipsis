process QUAST {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/reports/", pattern: "transposed_report.tsv", mode: "copy", saveAs: {"Quast_report.tsv"}

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

