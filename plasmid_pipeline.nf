assemblyfiles = Channel
	.fromPath(params.assemblies)
	.map { file -> tuple(file.baseName, file) }

readfiles = Channel
	.fromFilePairs(params.reads)



process MOBSUITE {
        publishDir "${params.outdir}/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.outdir}/mobsuite/mobtyper_reports", pattern: "*_mobtyper_plasmid_*.fasta_report.txt", mode: "copy"
 	publishDir "${params.outdir}/mobsuite", pattern: "*mob_recon.log", mode: "copy"

        input:
        set datasetID, file(datasetFile) from assemblyfiles
        
	output:
        file("*")
        file("*plasmid*fasta") into (resFasta, virFasta, plasFasta, prokkaFasta, mapFasta) mode flatten
	file("*_mobtyper_plasmid_*.fasta_report.txt") into (mobreportPmlst, mobreportNCBI)
        
	errorStrategy 'ignore'
        
        """
        mob_recon --infile ${datasetFile} -c --debug --run_typer --outdir . &> mob_recon.log
        rename '' ${datasetID}_ *
        """
}

process RESFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        publishDir "${params.outdir}/resfinder", pattern: "*results_tab.tsv", mode: "copy"
	publishDir "${params.outdir}/resfinder", pattern: "*resfinder.log", mode: "copy"

	input:
	tuple x, file("*plasmid*fasta") from resFasta
        
	output:
        file("*")
        
	"""
        python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i $x -o . -x -p /cluster/projects/nn9305k/src/resfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> resfinder.log
        rename '' "$x.baseName"_ *
	"""
}

process VIRFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        publishDir "${params.outdir}/virfinder", pattern: "*results_tab.tsv", mode: "copy"
	publishDir "${params.outdir}/virfinder", pattern: "*virfinder.log", mode: "copy"

        input:
        tuple x, file("*plasmid*fasta") from virFasta

        output:
        file("*")

        """
        python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i $x -o . -x -p /cluster/projects/nn9305k/src/virulencefinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> virfinder.log
        rename '' "$x.baseName"_ *
        """
}

process PLASMIDFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        publishDir "${params.outdir}/plasmidfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/plasmidfinder", pattern: "*plasmidfinder.log", mode: "copy"

        input:
        tuple x, file("*plasmid*fasta") from plasFasta

        output:
        file("*")

        """
        python /cluster/projects/nn9305k/src/plasmidfinder/plasmidfinder.py -i $x -o . -x -p /cluster/projects/nn9305k/src/plasmidfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> plasmidfinder.log
        rename '' "$x.baseName"_ *
        """
}
/*
process PROKKA {
	publishDir "${params.outdir}/prokka/gffs", pattern: "*.gff", mode: "copy"
	publishDir "${params.outdir}/prokka/gbk", pattern: "*.gbk", mode: "copy"
	publishDir "${params.outdir}/prokka/logs", pattern: "*.log", mode: "copy"

	input:
	tuple x, file("*plasmid*fasta") from prokkaFasta

	output:
	file("*")

	"""
	prokka --addgenes --compliant --prefix "$x.baseName" $x
	sed '/^##FASTA/Q' *.gff > "$x.baseName"_noseq.gff
	"""	
}
*/

process MAP {
	publishDir "${params.outdir}/mapping/mapped", pattern: "*txt", mode: "copy"
	tag { pair_id }

	input:
	set pair_id, file(rawreads) from readfiles
	tuple x, file("*plasmid*fasta") from mapFasta

	output:
	file("*")

	"""
	echo $x $rawreads > ${x}_results.txt	
	"""
}
