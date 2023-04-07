process BATCH_NORMALIZATION {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "conda-forge::r-base bioconda::bioconductor-deseq2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'quay.io/biocontainers/bioconductor-deseq2:1.34.0--r41h399db7b_0' }"

    input:
    path counts
    path samplesheet

    output:
    path "normalized.tsv", emit: normalized

    script:
    """
    cp $counts normalized.tsv
    """
}