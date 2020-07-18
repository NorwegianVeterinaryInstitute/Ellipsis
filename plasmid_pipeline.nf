log.info "================================================="
log.info " Plasmid Reconstruction and Improvement Pipeline"
log.info "================================================="
log.info "Assemblies		: ${params.assemblies}"
log.info "Reads                 : ${params.reads}"
log.info "AMR database		: ${params.amr_db}"
log.info "Virulence database	: ${params.vir_db}"
log.info "Output directory	: ${params.outdir}"
log.info "================================================="
log.info ""

// Channels
Channel
        .fromPath(params.assemblies)
        .map { file -> tuple(file.baseName, file) }
        .set { assemblyfiles }

Channel
        .fromFilePairs(params.reads, flat: true)
        .set { readfiles }

Channel
	.fromPath(params.amr_db)
	.set { amrDB_ch }

Channel
	.fromPath(params.vir_db)
	.set { virDB_ch }


// Pipeline
process MOBSUITE {
        publishDir "${params.outdir}/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.outdir}/mobsuite/mobtyper_reports", pattern: "*_mobtyper_plasmid_*.fasta_report.txt", mode: "copy"
 	publishDir "${params.outdir}/mobsuite", pattern: "*mob_recon.log", mode: "copy"
	
        tag "$datasetID"

        input:
        set datasetID, file(datasetFile) from assemblyfiles

        output:
        file("*")
        tuple datasetID, file("*mobtyper_plasmid*report.txt") into mobreport
        tuple datasetID, file("*plasmid*fasta") into plasmidFasta

        errorStrategy 'ignore'

        """
        mob_recon --infile ${datasetFile} -c --debug --run_typer --outdir . &> mob_recon.log
        rename '' "$datasetID"_ *
        """
}
plasmidFasta.transpose(remainder: true)
            .into{resFasta; virFasta; plasFasta; mapFasta}

process RESFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/resfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/resfinder", pattern: "*resfinder.log", mode: "copy"

	tag "$datasetID"

        input:
        set datasetID, file(plasmid) from resFasta

        output:
        file("*")

        """
        python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/resfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> resfinder.log
        rename '' "$plasmid.baseName"_ *
        """
}

process VIRFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/virfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/virfinder", pattern: "*virfinder.log", mode: "copy"

	tag "$datasetID"

        input:
        set datasetID, file(plasmid) from virFasta

        output:
        file("*")

        """
        python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/virulencefinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> virfinder.log
        rename '' "$plasmid.baseName"_ *
        """
}

process PLASFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"
        
	publishDir "${params.outdir}/plasmidfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.outdir}/plasmidfinder", pattern: "*plasmidfinder.log", mode: "copy"

	tag "$datasetID"	

        input:
        set datasetID, file(plasmid) from plasFasta

        output:
        file("*")

        """
        python /cluster/projects/nn9305k/src/plasmidfinder/plasmidfinder.py -i $plasmid -o . -x -p /cluster/projects/nn9305k/src/plasmidfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> plasmidfinder.log
        rename '' "$plasmid.baseName"_ *
        """
}

readfiles.combine(mapFasta, by: 0)
         .set {bbmapCh}

process BBDUK {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/BBTools"

        publishDir "${params.outdir}/bbduk", pattern: "*mapped*", mode: "copy"
        publishDir "${params.outdir}/bbduk", pattern: "*matched*", mode: "copy"
        publishDir "${params.outdir}/bbduk", pattern: "*log", mode: "copy"

	tag "$datasetID"

        input:
        tuple datasetID, file(R1), file(R2), file(plasmid) from bbmapCh

        output:
        file("*")
        tuple datasetID, file("*1_matched*"), file("*2_matched*") into (mappedReads_ch, aribaReads_ch)

        """
        bbduk.sh ref=$plasmid in=$R1 in2=$R2 out=${plasmid}_unmapped.fastq.gz outm1=${plasmid}_1_matched.fastq.gz outm2=${plasmid}_2_matched.fastq.gz &> ${plasmid}_bbduk.log
        """
}

process UNICYCLER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/unicycler"

        publishDir "${params.outdir}/unicycler/", pattern: "*assembly.fasta", mode: "copy"
        publishDir "${params.outdir}/unicycler/", pattern: "*unicycler.log", mode: "copy"

	tag "$datasetID"
	label 'heavy'

        input:
        tuple datasetID, file(R1), file(R2) from mappedReads_ch

        output:
        file("*")

	when:
	!params.complete

        """
        unicycler -1 $R1 -2 $R2 -o . --verbosity 2 --keep 2 --mode $params.mode --threads $task.cpus
        rename '' "$R1.baseName"_ *
        """
}

process ARIBA_RES {

	input:
	tuple datasetID, file(R1), file(R2) from aribaReads_ch
	file "db" from amrDB_ch

	"""
	echo ariba run --threads $task.cpus $db $R1 $R2 ${datasetID}_ariba
	"""
}
