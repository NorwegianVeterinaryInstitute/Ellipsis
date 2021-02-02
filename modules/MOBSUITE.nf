process MOB_RECON {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/mobsuite"

        publishDir "${params.out_dir}/results/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite", pattern: "*mob_recon.log", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/contig_reports", pattern: "*contig_report.txt", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/chromosome_fasta", pattern: "*chromosome.fasta", mode: "copy"

        tag "$datasetID"

        input:
        tuple val(datasetID), file(assembly)

        output:
        file("*")
        tuple val(datasetID), file("*plasmid*fasta"), emit: plasmidFasta
        tuple val(datasetID), file("*chromosome.fasta"), emit: chromFasta
        path "*contig_report.txt", emit: R_cont
	path "*mobtyper_results.txt", emit: R_mob

        errorStrategy 'ignore'

        script:
        """
        mob_recon -i $assembly -u -t --debug -n 4 -o results &> mob_recon.log
	mv results/* .
	rename '' "$datasetID"_ *
        """
}
