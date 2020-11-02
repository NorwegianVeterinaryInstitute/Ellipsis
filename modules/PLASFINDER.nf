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

