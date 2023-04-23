process DESEQ2 {
    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path counts
    path samplesheet
    tuple val(group1), val(group2)

    output:
    tuple val(group1), val(group2), path("*.tsv")

    script:
    """
    deseq.R --metadata $samplesheet --input $counts --group1 $group1 --group2 $group2
    """
}

process DESEQ2_GROUP {
    conda "pandas"
    container "tfprio-python"

    input:
        path(files)

    output:
        path("deseq_combined.tsv")

    script:
        """
        deseq2_group.py -d $files -o deseq_combined.tsv
        """
}