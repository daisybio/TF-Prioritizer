process STARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::stare-abc"
    container "biocontainers/stare-abc:1.0.4--haf6292c_1"
    
    input:
        tuple val(meta), path(candidate_regions)
        path(psem)
        path(annotation)
        path(genome_fasta)
        val(window_size)
        path(blacklist)
        val(decay)

    output:
        tuple val(meta), path("${meta.id}/Gene_TF_matrices/${meta.id}_TF_Gene_Affinities.txt"), emit: affinities

    script:
    """
        STARE.sh -c ${task.cpus} -a ${annotation} -g ${genome_fasta} -p ${psem} -b ${candidate_regions} -w ${window_size} -x ${blacklist} -e ${decay} -o ${meta.id}
        gzip -fd ${meta.id}/Gene_TF_matrices/${meta.id}_TF_Gene_Affinities.txt.gz
    """
}