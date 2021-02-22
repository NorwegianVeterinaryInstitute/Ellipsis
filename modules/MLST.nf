process MLST_READS {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/stringmlst"

        publishDir "${params.out_dir}/results/stringmlst", pattern: "*stringmlst_report.tsv", mode: "copy"

        input:
        tuple val(datasetID), file(R1), file(R2)

        output:
        file("*")
        path "*stringmlst_report.tsv", emit: R_stringmlst

        script:
        """
        stringMLST.py --predict -1 $R1 -2 $R2 -p --prefix ${params.stringmlst_path}/${params.stringmlst_db} -k $params.kmersize -o ${datasetID}_stringmlst_report.tsv
        """
}

process MLST_ASSEMBLY {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/mlst"

        publishDir "${params.out_dir}/results/ariba", pattern: "*assembly_mlst_report.tsv", mode: "copy"

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

