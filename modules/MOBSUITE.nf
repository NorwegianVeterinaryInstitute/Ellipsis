process MOB_RECON {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Mobsuite"

        publishDir "${params.out_dir}/results/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/mobtyper_reports", pattern: "*_mobtyper_plasmid_*.fasta_report.txt", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite", pattern: "*mob_recon.log", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/repetitive_blast_reports", pattern: "*repetitive_blast_report.txt", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/contig_reports", pattern: "*contig_report.txt", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/chromosome_fasta", pattern: "*chromosome.fasta", mode: "copy"

        tag "$datasetID"

        input:
        tuple val(datasetID), file(assembly)

        output:
        file("*")
        tuple val(datasetID), file("*mobtyper_plasmid*report.txt"), emit: mobreport
        tuple val(datasetID), file("*plasmid*fasta"), emit: plasmidFasta
        tuple val(datasetID), file("*chromosome.fasta"), emit: chromFasta
        path "*mobtyper_plasmid*report.txt", emit: R_mob

        errorStrategy 'ignore'

        script:
        """
        mob_recon --infile $assembly -c --debug --run_typer --outdir . &> mob_recon.log
        rename '' "$datasetID"_ *
        """
}

