process MERGE_BINDING_BED {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(hm), val(group), val(tf), path(beds)

    output:
    tuple val(hm), val(group), val(tf), path("${hm}_${group}_${tf}_merged.bed")

    script:
    if (beds.getClass().isArray() && beds.size() > 1) {
        """
        bedtools merge -c 4 -o mean -i ${beds.join(' ')} > ${hm}_${group}_${tf}_merged.bed
        """
    } else {
        """
        ln -s ${beds[0]} ${hm}_${group}_${tf}_merged.bed
        """
    }
}