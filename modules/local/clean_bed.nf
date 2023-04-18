process CLEAN_BED {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
        tuple val(sample), path(bed)

    output:
        tuple val(sample), path("${sample}_sorted.bed")

    script:
        """
        cut -f1-3 ${bed} | bedtools sort | awk -v sample=${sample} '{print \$0"\t"sample}' > ${sample}_sorted.bed
        """
}