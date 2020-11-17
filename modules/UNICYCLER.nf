process UNICYCLER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/unicycler"

        publishDir "${params.out_dir}/results/unicycler/", pattern: "*assembly.fasta", mode: "copy"
        publishDir "${params.out_dir}/results/unicycler/", pattern: "*unicycler.log", mode: "copy"

        tag "$datasetID"
        label 'longtime'

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
        tuple val(datasetID), path {"*_assembly.fasta"}, emit: new_assemblies
        path "*_assembly.fasta", emit: quast_ch

        """
        unicycler -1 $R1 -2 $R2 -o . --verbosity 2 --keep 2 --mode $params.mode --threads $task.cpus --min_fasta_length $params.min_fasta_length
        rename '' "$datasetID"_ *
        """
}

process UNICYCLER_HYBRID {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/unicycler"

        publishDir "${params.out_dir}/results/unicycler/", pattern: "*assembly.fasta", mode: "copy"
        publishDir "${params.out_dir}/results/unicycler/", pattern: "*unicycler.log", mode: "copy"

        tag "$datasetID"
        label 'hybrid'

        input:
        tuple val(datasetID), file(R1), file(R2), file(longreads)

        output:
        file("*")
        tuple val(datasetID), path {"*_assembly.fasta"}, emit: new_assemblies
        path "*_assembly.fasta", emit: quast_ch

        """
        unicycler -1 $R1 -2 $R2 -l $longreads -o . --verbosity 2 --keep 2 --mode $params.mode --threads $task.cpus --min_fasta_length $params.min_fasta_length
        rename '' "$datasetID"_ *
        """
}

