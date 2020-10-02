log.info "=================================================="
log.info "=================================================="
log.info "                                                  "
log.info "****                 ELLIPSIS                 ****"
log.info "                                                  "
log.info "****            Plasmid annotation            ****"
log.info "****         and improvement pipeline         ****"
log.info "                                                  "
log.info "=================================================="
log.info "=================================================="

nextflow.enable.dsl=2

// Processes

process MOB_RECON {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Mobsuite"

        publishDir "${params.out_dir}/mobsuite/plasmid_fasta", pattern: "*plasmid*fasta", mode: "copy"
        publishDir "${params.out_dir}/mobsuite/mobtyper_reports", pattern: "*_mobtyper_plasmid_*.fasta_report.txt", mode: "copy"
        publishDir "${params.out_dir}/mobsuite", pattern: "*mob_recon.log", mode: "copy"
        publishDir "${params.out_dir}/mobsuite/repetitive_blast_reports", pattern: "*repetitive_blast_report.txt", mode: "copy"
        publishDir "${params.out_dir}/mobsuite/contig_reports", pattern: "*contig_report.txt", mode: "copy"
        publishDir "${params.out_dir}/mobsuite/chromosome_fasta", pattern: "*chromosome.fasta", mode: "copy"

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

process RESFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"

        publishDir "${params.out_dir}/resfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.out_dir}/resfinder", pattern: "*resfinder.log", mode: "copy"

        tag "$fasta.baseName"

        input:
        tuple val(datasetID), file(fasta)

        output:
        file("*")
        path "*results_tab.tsv", emit: R_res

        """
        python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i $fasta -o . -x -p $params.resfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> resfinder.log
        rename '' "$fasta.baseName"_"resfinder"_ *
        """
}

process VIRFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"

        publishDir "${params.out_dir}/virfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.out_dir}/virfinder", pattern: "*virfinder.log", mode: "copy"

        tag "$fasta.baseName"

        input:
        tuple val(datasetID), file(fasta)

        output:
        file("*")
        path "*results_tab.tsv", emit: R_vir

        """
        python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i $fasta -o . -x -p $params.virfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> virfinder.log
        rename '' "$fasta.baseName"_"virfinder"_ *
        """
}


process PLASFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"

        publishDir "${params.out_dir}/plasmidfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.out_dir}/plasmidfinder", pattern: "*plasmidfinder.log", mode: "copy"

        tag "$fasta.baseName"

        input:
        tuple val(datasetID), file(fasta)

        output:
        file("*")
        path "*results_tab.tsv", emit: R_plas

        """
        python /cluster/projects/nn9305k/src/plasmidfinder/plasmidfinder.py -i $fasta -o . -x -p $params.plasfinder_db -mp /cluster/software/BLAST+/2.8.1-foss-2018b/bin/blastn &> plasmidfinder.log
        rename '' "$fasta.baseName"_"plasfinder"_ *
        """
}


process PROKKA {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/prokka/${fasta.baseName}", pattern: "*", mode: "copy"

        tag "$fasta.baseName"

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


process ARIBA {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"
	
	publishDir "${params.out_dir}/ariba", pattern: "*ariba_report.tsv", mode: "copy"
	publishDir "${params.out_dir}/ariba", pattern: "*_ariba.log", mode: "copy"

	input:
	tuple val(datasetID), file(R1), file(R2)
	path amrdb

	output:
	file("*")
	path "*ariba_report.tsv", emit: R_ariba

	script:
	"""
	ariba run --threads $task.cpus $amrdb $R1 $R2 results &> ${datasetID}_ariba.log
	cp results/report.tsv ${datasetID}_ariba_report.tsv
	"""
}


process TRIM {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/Trimgalore"

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
	tuple val(datasetID), path {"*val_1.fq.gz"}, path {"*val_2.fq.gz"}, emit: trim_reads 

	script:
        """
	trim_galore -o . --paired -trim1 $R1 $R2 &> ${datasetID}_trimgalore.log
        """
}


process UNICYCLER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/unicycler"

        publishDir "${params.out_dir}/unicycler/", pattern: "*assembly.fasta", mode: "copy"
        publishDir "${params.out_dir}/unicycler/", pattern: "*unicycler.log", mode: "copy"

        tag "$datasetID"
        label 'heavy'

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
	tuple val(datasetID), path {"*_assembly.fasta"}, emit: new_assemblies
	path "*_assembly.fasta", emit: quast_ch

        """
        unicycler -1 $R1 -2 $R2 -o . --verbosity 2 --keep 2 --mode $params.mode --threads $task.cpus
        rename '' "$datasetID"_ *
        """
}


process QUAST {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

	publishDir "${params.out_dir}/reports/", pattern: "transposed_report.tsv", mode: "copy", saveAs: {"quast_report.tsv"}

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

process REPORT {
        module 'R/4.0.0-foss-2020a'

        publishDir "${params.out_dir}/reports", pattern: "total_report.txt", mode: "copy"
	publishDir "${params.out_dir}/reports", pattern: "summary_report.txt", mode: "copy"
	publishDir "${params.out_dir}/reports", pattern: "resfinder_report.txt", mode: "copy"
	publishDir "${params.out_dir}/reports", pattern: "virfinder_report.txt", mode: "copy"
	publishDir "${params.out_dir}/reports", pattern: "plasmidfinder_report.txt", mode: "copy"

        input:
        file("*")

        output:
        file("*")

        """
        Rscript $baseDir/bin/collate_data.R
        """
}

// Workflows

workflow ELLIPSIS_ASSEMBLY {
	Channel
                .fromFilePairs(params.reads, flat: true,  checkIfExists: true)
		.set { readfiles }

	aribadb = Channel
		.value(params.ariba_db)

	ARIBA(readfiles, aribadb)
	TRIM(readfiles)
	UNICYCLER(TRIM.out.trim_reads)
	QUAST(UNICYCLER.out.quast_ch.collect())
	MOB_RECON(UNICYCLER.out.new_assemblies)

	if (!params.chrom) {
		MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                	.set { fasta_ch }
	}
	if (params.chrom) {
		MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                	.mix(MOB_RECON.out.chromFasta.transpose(remainder: true))
                	.set { fasta_ch }
	}

	RESFINDER(fasta_ch)
        VIRFINDER(fasta_ch)
        PLASFINDER(fasta_ch)
        PROKKA(fasta_ch)

	RESFINDER.out.R_res.collect()
                .mix(VIRFINDER.out.R_vir)
                .mix(PLASFINDER.out.R_plas)
                .mix(PROKKA.out.R_prokka)
                .mix(MOB_RECON.out.R_mob)
		.mix(ARIBA.out.R_ariba)
		.mix(QUAST.out.R_quast)
                .collect()
                .set { report_ch }

        REPORT(report_ch)
}

workflow ELLIPSIS_ANNOTATE {
	Channel
                .fromPath(params.assemblies, checkIfExists: true)
                .map { file -> tuple(file.baseName, file) }
                .set { assembly_ch }

	MOB_RECON(assembly_ch)

        if (!params.chrom) {
                MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                        .set { fasta_ch }
        }
        if (params.chrom) {
                MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                        .mix(MOB_RECON.out.chromFasta.transpose(remainder: true))
                        .set { fasta_ch }
        }

        RESFINDER(fasta_ch)
        VIRFINDER(fasta_ch)
        PLASFINDER(fasta_ch)
        PROKKA(fasta_ch)

/*
        RESFINDER.out.R_res.collect()
                .mix(VIRFINDER.out.R_vir)
                .mix(PLASFINDER.out.R_plas)
                .mix(PROKKA.out.R_prokka)
                .mix(MOB_RECON.out.R_mob)
                .collect()
                .set { report_ch }

        REPORT(report_ch) */
}



workflow {
if (params.assemble) {
	ELLIPSIS_ASSEMBLY()
	}

if (!params.assemble) {
	ELLIPSIS_ANNOTATE()
	}
}
