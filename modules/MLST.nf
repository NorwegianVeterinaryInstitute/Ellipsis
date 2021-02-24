process MLST {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/mlst"

        publishDir "${params.out_dir}/results/mlst", pattern: "*assembly_mlst_report.tsv", mode: "copy"

        input:
        tuple val(datasetID), file(assembly)

        output:
        file("*")
        path "*assembly_mlst_report.tsv", emit: R_mlst

        script:
        """
	mlst --legacy --scheme $params.mlst_species $assembly > ${datasetID}_assembly_mlst_report.tsv
        """
}

