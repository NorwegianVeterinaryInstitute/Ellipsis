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

