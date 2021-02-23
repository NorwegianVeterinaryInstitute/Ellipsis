process ARIBA_RES {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/ariba", pattern: "*ariba_resfinder_report.tsv", mode: "copy"
        publishDir "${params.out_dir}/results/ariba", pattern: "*_ariba_res.log", mode: "copy"

        input:
        tuple val(datasetID), file(R1), file(R2)
        path amrdb

        output:
        file("*")
        path "*ariba_resfinder_report.tsv", emit: R_aribares

        script:
        """
        ariba run --threads $task.cpus $amrdb $R1 $R2 results &> ${datasetID}_ariba_res.log
        cp results/report.tsv ${datasetID}_ariba_resfinder_report.tsv
        """
}

process ARIBA_VIR {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/ariba", pattern: "*ariba_virulence_report.tsv", mode: "copy"
        publishDir "${params.out_dir}/results/ariba", pattern: "*_ariba_vir.log", mode: "copy"

        input:
        tuple val(datasetID), file(R1), file(R2)
        path virdb

        output:
        file("*")
        path "*ariba_virulence_report.tsv", emit: R_aribavir

        script:
        """
        ariba run --threads $task.cpus $virdb $R1 $R2 results &> ${datasetID}_ariba_vir.log
        cp results/report.tsv ${datasetID}_ariba_virulence_report.tsv
        """
}

process ARIBA_MLST {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/bifrost"

        publishDir "${params.out_dir}/results/ariba", pattern: "*ariba_mlst_report.tsv", mode: "copy"
	publishDir "${params.out_dir}/results/ariba", pattern: "*ariba_mlst_detailed_report.tsv", mode: "copy"
        publishDir "${params.out_dir}/results/ariba", pattern: "*_ariba_mlst.log", mode: "copy"

        input:
        tuple val(datasetID), file(R1), file(R2)
        path mlstdb

        output:
        file("*")
        path "*ariba_mlst_report.tsv", emit: R_aribamlst

        script:
        """
        ariba run --threads $task.cpus $mlstdb $R1 $R2 results &> ${datasetID}_ariba_mlst.log
        cp results/report.tsv ${datasetID}_ariba_mlst_detailed_report.tsv
	cp results/mlst_report.tsv ${datasetID}_ariba_mlst_report.tsv
        """
}
