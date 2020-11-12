process RESFINDER {
        conda "/cluster/projects/nn9305k/src/miniconda/envs/cge_addons"

        publishDir "${params.out_dir}/results/resfinder", pattern: "*results_tab.tsv", mode: "copy"
        publishDir "${params.out_dir}/results/resfinder", pattern: "*resfinder.log", mode: "copy"

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

