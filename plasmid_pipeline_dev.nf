log.info "================================================="
log.info " Plasmid Reconstruction and Improvement Pipeline"
log.info "================================================="
log.info "Assemblies		: ${params.assemblies}"
log.info "Reads			: ${params.reads}"
log.info "AMR database		: ${params.amr_db}"
log.info "Virulence database	: ${params.vir_db}"
log.info "Output directory	: ${params.outdir}"
log.info "================================================="
log.info ""

// Channels
Channel
        .fromPath(params.assemblies, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .set { assemblyfiles }

Channel
        .fromFilePairs(params.reads, flat: true, checkIfExists: true)
        .set { readfiles }

Channel
	.value(params.amr_db)
	.set { amrDB_ch }

Channel
	.value(params.vir_db)
	.set { virDB_ch }


// Pipeline

process MOB_RECON {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/Mobsuite"

        publishDir "${params.outdir}/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.outdir}/mobsuite/mobtyper_reports", pattern: "*_mobtyper_plasmid_*.fasta_report.txt", mode: "copy"
 	publishDir "${params.outdir}/mobsuite", pattern: "*mob_recon.log", mode: "copy"
	publishDir "${params.outdir}/mobsuite/repetitive_blast_reports", pattern: "*repetitive_blast_report.txt", mode: "copy"
	publishDir "${params.outdir}/mobsuite/contig_reports", pattern: "*contig_report.txt", mode: "copy"
	publishDir "${params.outdir}/mobsuite/chromosome_fasta", pattern: "*chromosome.fasta", mode: "copy"	

        tag "$datasetID"

        input:
        tuple datasetID, file(datasetFile) from assemblyfiles

        output:
        file("*")
        tuple datasetID, file("*mobtyper_plasmid*report.txt") into mobreport
        tuple datasetID, file("*plasmid*fasta") into plasmidFasta
	file("*mobtyper_plasmid*report.txt") into R_mob

        errorStrategy 'ignore'

	script:
        """
        mob_recon --infile $datasetFile -c --debug --run_typer --outdir . &> mob_recon.log
	rename '' "$datasetID"_ *
        """
}

mobreport.transpose(remainder: true)
	 .set { mobreport_ch }

process MASH_ACCESSION {
	
	input:
	tuple datasetID, file(report) from mobreport_ch

	output:
	tuple env(plasmid), stdout into mash_acc

	"""
	get_accession.bash $report
	plasmid_name0=${report/_mobtyper}
	plasmid=${plasmid_name0%%.fasta_report.txt}
	"""
}

mash_acc.view()

/*
mash_acc.map { file -> tuple(file.baseName, file) }
	.set { mash_acc_ch }

process NCBIMAP {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/ncbiacc"

	input:
	tuple plasmid, file(acc) from mash_acc_ch

	output:
	file("*")
	
	script:
	mash_acc = acc.name
	"""
	ncbi-acc-download -m nucleotide -F fasta -p $plasmid $mash_acc
	"""
}

plasmidFasta.transpose(remainder: true)
            .into{resFasta; virFasta; plasFasta; mapFasta; prokkaFasta}

process RESFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/resfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/resfinder", pattern: "*resfinder.log", mode: "copy"

	tag "$plasmid.baseName"

        input:
        set datasetID, file(plasmid) from resFasta

        output:
        file("*")
	file("*results_tab.tsv") into R_res

        """
        python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/resfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> resfinder.log
        rename '' "$plasmid.baseName"_"resfinder"_ *
        """
}

process VIRFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/virfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/virfinder", pattern: "*virfinder.log", mode: "copy"

	tag "$plasmid.baseName"

        input:
        set datasetID, file(plasmid) from virFasta

        output:
        file("*")
	file("*results_tab.tsv") into R_vir

        """
        python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/virulencefinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> virfinder.log
        rename '' "$plasmid.baseName"_"virfinder"_ *
        """
}

process PLASFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/plasmidfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/plasmidfinder", pattern: "*plasmidfinder.log", mode: "copy"

	tag "$plasmid.baseName"	

        input:
        set datasetID, file(plasmid) from plasFasta

        output:
        file("*")
	file("*results_tab.tsv") into R_plas

        """
        python /cluster/projects/nn9305k/src/plasmidfinder/plasmidfinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/plasmidfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> plasmidfinder.log
        rename '' "$plasmid.baseName"_"plasfinder"_ *
        """
}

process PROKKA {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

	publishDir "${params.outdir}/prokka/${plasmid.baseName}", pattern: "*", mode: "copy"

	tag "$plasmid.baseName"

	input:
	set datasetID, file(plasmid) from prokkaFasta

	output:
	file("*")
	file("*.txt") into R_prokka

	"""
	prokka --addgenes --compliant --force --prefix $plasmid.baseName --outdir . $plasmid
	rename '' "prokka_report_" *
	"""
}

readfiles.combine(mapFasta, by: 0)
         .set { bbdukCh }

process BBDUK {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/BBTools"

        publishDir "${params.outdir}/bbduk", pattern: "*mapped*", mode: "copy"
        publishDir "${params.outdir}/bbduk", pattern: "*matched*", mode: "copy"
        publishDir "${params.outdir}/bbduk", pattern: "*log", mode: "copy"

	tag "$plasmid.baseName"
	label 'medium'

        input:
        tuple datasetID, file(R1), file(R2), file(plasmid) from bbdukCh

        output:
        file("*")
        tuple file(plasmid), file("*1_matched*"), file("*2_matched*") into (mappedReads_ch, aribaReads_ch)

	when:
        !params.complete

        """
        bbduk.sh ref=$plasmid in=$R1 in2=$R2 out=${plasmid.baseName}_unmapped.fastq.gz outm1=${plasmid.baseName}_1_matched.fastq.gz outm2=${plasmid.baseName}_2_matched.fastq.gz &> ${plasmid.baseName}_bbduk.log
        """
}

process UNICYCLER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/unicycler"

	publishDir "${params.outdir}/unicycler/", pattern: "*assembly.fasta", mode: "copy"
        publishDir "${params.outdir}/unicycler/", pattern: "*unicycler.log", mode: "copy"

	tag "$plasmid.baseName"
	label 'heavy'

        input:
        tuple file(plasmid), file(R1), file(R2) from mappedReads_ch

        output:
        file("*")

	when:
	!params.complete

        """
        unicycler -1 $R1 -2 $R2 -o . --verbosity 2 --keep 2 --mode $params.mode --threads $task.cpus
        rename '' "$plasmid.baseName"_ *
        """
}


process ARIBA_AMR {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

	publishDir "${params.outdir}/ariba/amr", pattern: "*report.tsv", mode: "copy"
	publishDir "${params.outdir}/ariba/amr", pattern: "*ariba.log", mode: "copy"

	tag "$plasmid.baseName"

	input:
	tuple plasmid, file(R1), file(R2) from aribaReads_ch
	val "db" from amrDB_ch

	output:
	file("*")
	tuple $plasmid.baseName, file("*_report.tsv") into collateAribaAmr

	"""
	ariba run --threads $task.cpus --assembly_cov 200 --assembler spades --spades_mode wgs --verbose $db $R1 $R2 ${plasmid.baseName}_ariba &> ariba.log
	cp ${plasmid.baseName}_ariba/report.tsv .
	rename '' "$plasmid.baseName"_ *
	"""
}


process COLLATE {
	module 'R/4.0.0-foss-2020a'

	publishDir "${params.outdir}", pattern: "total_report.txt", mode: "copy"

	input:
	file("*") from R_mob.collect()
	file("*") from R_res.collect()
	file("*") from R_vir.collect()
	file("*") from R_plas.collect()
	file("*") from R_prokka.collect()

	output:
	file("*")

	"""
	Rscript $baseDir/bin/collate_data.R
	"""
}

*/
