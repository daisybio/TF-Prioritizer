process MIX_SAMPLES {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
        val(min_occurrence)
        tuple val(hm), val(group), path(peak_files)

    output:
        tuple val(hm), val(group), val("${hm}_${group}_merged"), path("${hm}_${group}.bed")

    script:
        """
        cat ${peak_files.join(' ')} | bedtools sort | bedtools merge -c 4 -o count_distinct | awk -v cutoff=$min_occurrence '\$NF >= cutoff' | cut -f1-3 | bedtools sort > ${hm}_${group}.bed
        """
}