process PROKKA {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/prokka/${fasta.baseName}", pattern: "*", mode: "copy"

        tag "$fasta.baseName"

	label 'longtime'

        input:
        tuple val(datasetID), file(fasta)

        output:
        file("*")
        path "*.txt", emit: R_prokka

        """
        prokka --addgenes --compliant --force --cpus $task.cpus --prefix $fasta.baseName --outdir . $params.prokka_additional $fasta
        rename '' "prokka_report_" *
        """
}

