process MOB_RECON {
	conda (params.enable_conda ? 'bioconda::mob_suite=3.1.8' : null)
        container 'quay.io/biocontainers/mob_suite:3.1.8--pyhdfd78af_1'

        publishDir "${params.out_dir}/results/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite", pattern: "*mob_recon.log", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/contig_reports", pattern: "*contig_report.txt", mode: "copy"
        publishDir "${params.out_dir}/results/mobsuite/chromosome_fasta", pattern: "*chromosome.fasta", mode: "copy"

        input:
        tuple val(datasetID), file(assembly)

        output:
        file("*")
        tuple val(datasetID), file("*plasmid*fasta"), optional: true, emit: plasmidFasta
        tuple val(datasetID), file("*chromosome.fasta"), emit: chromFasta
        path "*contig_report.txt", emit: R_cont
	path "*mobtyper_results.txt", optional: true, emit: R_mob

        script:
        """
        mob_recon -i $assembly -u -t --debug -n 4 -o results &> mob_recon.log
	mv results/* .
	rename '' "$datasetID"_ *
        """
}
