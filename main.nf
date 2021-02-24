
log.info "=================================================="
log.info "=================================================="
log.info "                                                  "
log.info "****                 ELLIPSIS                 ****"
log.info "                                                  "
log.info "****           Plasmid identification,        ****"
log.info "****       annotation and typing pipeline     ****"
log.info "                                                  "
log.info "=================================================="
log.info "=================================================="

nextflow.enable.dsl=2

// Workflows
workflow ELLIPSIS_HYBRID {
	// Set channels
	Channel
                .fromFilePairs(params.reads, flat: true,  checkIfExists: true)
                .set { readfiles_ch }

	Channel
		.fromPath(params.longreads, checkIfExists: true)
		.map { file -> tuple(file.simpleName, file) }
		.set { longreads_ch }

	Channel
                .value(params.ariba_resdb)
		.set { aribaresdb }

        Channel
                .value(params.ariba_virdb)
		.set { aribavirdb }

	Channel
                .value(params.ariba_mlstdb)
                .set { aribamlstdb }

	Channel
		.value(params.sequencer)
		.set { seq_type }

	Channel
		.value(params.illumina_filtering)
		.set { filt_illu }

	// Read QC
	FASTQC(readfiles_ch)
	NANOPLOT(longreads_ch)
	MULTIQC_PRE(FASTQC.out.fastqc_reports.collect())

	// Trimming and filtering
	if (params.trim) {
		TRIM(readfiles_ch)
		FASTQC_POST(TRIM.out.trim_reads)
		MULTIQC_POST(FASTQC_POST.out.fastqc_reports.collect())
		CANU(seq_type, longreads_ch)

		TRIM.out.trim_reads
                	.join(CANU.out.canu_output, by: 0)
                        .set { pre_filt_reads_ch }		

		FILTLONG(filt_illu, pre_filt_reads_ch)

		TRIM.out.trim_reads
                       	.join(FILTLONG.out.filtered_longreads, by: 0)
                       	.set { all_reads_ch }
	}
	if (!params.trim) {
		readfiles_ch
			.join(longreads_ch)
			.set { pre_filt_reads_ch }

		FILTLONG(filt_illu, pre_filt_reads_ch)
	
		readfiles_ch
			.join(FILTLONG.out.filtered_longreads, by: 0)
			.set { all_reads_ch }
	}

	UNICYCLER_HYBRID(all_reads_ch)
	QUAST(UNICYCLER_HYBRID.out.quast_ch.collect())        

	// Plasmid analyses
	MOB_RECON(UNICYCLER_HYBRID.out.new_assemblies)

        if (!params.chrom) {
                MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                        .set { fasta_ch }
        }
        if (params.chrom) {
                MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                        .mix(MOB_RECON.out.chromFasta.transpose(remainder: true))
                        .set { fasta_ch }
        }

        ARIBA_RES(readfiles_ch, aribaresdb)
        ARIBA_VIR(readfiles_ch, aribavirdb)
	ARIBA_MLST(readfiles_ch, aribamlstdb)
	MLST(UNICYCLER_HYBRID.out.new_assemblies)
        RESFINDER(fasta_ch)
        VIRFINDER(fasta_ch)
        PLASFINDER(fasta_ch)
        PROKKA(fasta_ch)

        RESFINDER.out.R_res.collect()
                .mix(VIRFINDER.out.R_vir)
                .mix(PLASFINDER.out.R_plas)
                .mix(PROKKA.out.R_prokka)
                .mix(MOB_RECON.out.R_mob)
		.mix(MOB_RECON.out.R_cont)
                .mix(ARIBA_RES.out.R_aribares)
                .mix(ARIBA_VIR.out.R_aribavir)
		.mix(ARIBA_MLST.out.R_aribamlst)
		.mix(MLST.out.R_mlst)
                .mix(QUAST.out.R_quast)
                .collect()
                .set { report_ch }

        Channel
                .value("true")
                .set { run_ariba_report }

        REPORT(report_ch, run_ariba_report)
}

workflow ELLIPSIS_ASSEMBLY {
	Channel
                .fromFilePairs(params.reads, flat: true,  checkIfExists: true)
		.set { readfiles_ch }

	Channel
		.value(params.ariba_resdb)
		.set { aribaresdb }

	Channel
                .value(params.ariba_virdb)
		.set { aribavirdb }

	Channel
                .value(params.ariba_mlstdb)
                .set { aribamlstdb }

	// Read QC
        FASTQC(readfiles_ch)
        MULTIQC_PRE(FASTQC.out.fastqc_reports.collect())

	// Assembly
	TRIM(readfiles_ch)
	FASTQC_POST(TRIM.out.trim_reads)
        MULTIQC_POST(FASTQC_POST.out.fastqc_reports.collect())
	UNICYCLER(TRIM.out.trim_reads)
	QUAST(UNICYCLER.out.quast_ch.collect())
	MOB_RECON(UNICYCLER.out.new_assemblies)

	// Plasmid analyses
	if (!params.chrom) {
		MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                	.set { fasta_ch }
	}
	if (params.chrom) {
		MOB_RECON.out.plasmidFasta.transpose(remainder: true)
                	.mix(MOB_RECON.out.chromFasta.transpose(remainder: true))
                	.set { fasta_ch }
	}

	ARIBA_RES(readfiles_ch, aribaresdb)
        ARIBA_VIR(readfiles_ch, aribavirdb)
	ARIBA_MLST(readfiles_ch, aribamlstdb)
	MLST(UNICYCLER.out.new_assemblies)
	RESFINDER(fasta_ch)
        VIRFINDER(fasta_ch)
        PLASFINDER(fasta_ch)
        PROKKA(fasta_ch)

	RESFINDER.out.R_res.collect()
                .mix(VIRFINDER.out.R_vir)
                .mix(PLASFINDER.out.R_plas)
                .mix(PROKKA.out.R_prokka)
                .mix(MOB_RECON.out.R_mob)
		.mix(MOB_RECON.out.R_cont)
		.mix(ARIBA_RES.out.R_aribares)
		.mix(ARIBA_VIR.out.R_aribavir)
		.mix(ARIBA_MLST.out.R_aribamlst)
		.mix(MLST.out.R_mlst)
		.mix(QUAST.out.R_quast)
                .collect()
                .set { report_ch }
	
	Channel
		.value("true")
		.set { run_ariba_report }

        REPORT(report_ch, run_ariba_report)
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

	MLST(assembly_ch)
        RESFINDER(fasta_ch)
        VIRFINDER(fasta_ch)
        PLASFINDER(fasta_ch)
        PROKKA(fasta_ch)

        RESFINDER.out.R_res.collect()
                .mix(VIRFINDER.out.R_vir)
                .mix(PLASFINDER.out.R_plas)
                .mix(PROKKA.out.R_prokka)
                .mix(MOB_RECON.out.R_mob)
		.mix(MOB_RECON.out.R_cont)
		.mix(MLST.out.R_mlst)
                .collect()
                .set { report_ch }

	Channel
                .value("false")
                .set { run_ariba_report }

        REPORT(report_ch, run_ariba_report)
}



workflow {
if (params.track == "hybrid") {
	include { FASTQC; FASTQC as FASTQC_POST } from "${params.module_dir}/FASTQC.nf"
	include { MULTIQC_PRE; MULTIQC_POST } from "${params.module_dir}/MULTIQC.nf"
	include { NANOPLOT } from "${params.module_dir}/NANOPLOT.nf"
	include { TRIM } from "${params.module_dir}/TRIM.nf"
	include { CANU } from "${params.module_dir}/CANU.nf"
	include { FILTLONG } from "${params.module_dir}/FILTLONG.nf"
	include { UNICYCLER_HYBRID } from "${params.module_dir}/UNICYCLER.nf"
        include { QUAST } from "${params.module_dir}/QUAST.nf"
        include { MOB_RECON } from "${params.module_dir}/MOBSUITE.nf"
	include { ARIBA_RES;ARIBA_VIR;ARIBA_MLST } from "${params.module_dir}/ARIBA.nf"
	include { MLST } from "${params.module_dir}/MLST.nf"
	include { RESFINDER } from "${params.module_dir}/RESFINDER.nf"
        include { VIRFINDER } from "${params.module_dir}/VIRFINDER.nf"
        include { PLASFINDER } from "${params.module_dir}/PLASFINDER.nf"
        include { PROKKA } from "${params.module_dir}/PROKKA.nf"
        include { REPORT } from "${params.module_dir}/REPORT.nf"

	ELLIPSIS_HYBRID()
}

if (params.track == "short_assembly") {
	include { FASTQC; FASTQC as FASTQC_POST } from "${params.module_dir}/FASTQC.nf"
        include { MULTIQC_PRE; MULTIQC_POST } from "${params.module_dir}/MULTIQC.nf"
	include { ARIBA_RES;ARIBA_VIR;ARIBA_MLST } from "${params.module_dir}/ARIBA.nf"
	include { MLST } from "${params.module_dir}/MLST.nf"
	include { TRIM } from "${params.module_dir}/TRIM.nf"
	include { UNICYCLER } from "${params.module_dir}/UNICYCLER.nf"
	include { QUAST } from "${params.module_dir}/QUAST.nf"
	include { MOB_RECON } from "${params.module_dir}/MOBSUITE.nf"
	include { RESFINDER } from "${params.module_dir}/RESFINDER.nf"
	include { VIRFINDER } from "${params.module_dir}/VIRFINDER.nf"
	include { PLASFINDER } from "${params.module_dir}/PLASFINDER.nf"
	include { PROKKA } from "${params.module_dir}/PROKKA.nf"
	include { REPORT } from "${params.module_dir}/REPORT.nf"

	ELLIPSIS_ASSEMBLY()
	}

if (params.track == "no_assembly") {
        include { MOB_RECON } from "${params.module_dir}/MOBSUITE.nf"
        include { RESFINDER } from "${params.module_dir}/RESFINDER.nf"
        include { VIRFINDER } from "${params.module_dir}/VIRFINDER.nf"
        include { PLASFINDER } from "${params.module_dir}/PLASFINDER.nf"
	include { MLST } from "${params.module_dir}/MLST.nf"
        include { PROKKA } from "${params.module_dir}/PROKKA.nf"
        include { REPORT } from "${params.module_dir}/REPORT.nf"

	ELLIPSIS_ANNOTATE()
	}
}

workflow.onComplete {
	log.info "".center(60, "=")
	log.info "Ellipsis Complete!".center(60)
	log.info "Output directory: $params.out_dir".center(60)
	log.info "Duration: $workflow.duration".center(60)
	log.info "Command line: $workflow.commandLine".center(60)
	log.info "Nextflow version: $workflow.nextflow.version".center(60)
	log.info "".center(60, "=")
}

workflow.onError {
	println "Pipeline execution stopped with the following message: ${workflow.errorMessage}".center(60, "=")
}
