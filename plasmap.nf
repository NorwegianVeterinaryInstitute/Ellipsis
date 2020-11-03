log.info "=================================================="
log.info "                                                  "
log.info "****         Plasmid Mapping Pipeline         ****"
log.info "                                                  "
log.info "=================================================="

nextflow.enable.dsl=2

// Workflow
workflow PLASMAP {
	Channel
		.fromPath(params.ref)
		.map { file -> tuple(file.baseName, file) }
		.set { refs_ch }

	Channel
                .fromFilePairs(params.reads, flat: true,  checkIfExists: true)
                .combine(refs_ch)
		.set { bwa_ch }

	BWA_MEM_SAM(bwa_ch)
	COV_CALC(BWA_MEM_SAM.out.bam_ch)

	COV_CALC.out.R_ch.collect()
		.set { R_reports_ch }

	REPORT_PLASMAP(R_reports_ch)
}

workflow {
	include { BWA_MEM_SAM } from "${params.module_dir}/BWA.nf"
	include { COV_CALC } from "${params.module_dir}/COVCALC.nf"
	include { REPORT_PLASMAP } from "${params.module_dir}/REPORT.nf"

	PLASMAP()
}
