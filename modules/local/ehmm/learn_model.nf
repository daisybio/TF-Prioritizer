process LEARN_MODEL {
    tag "$meta.id"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        tuple val(meta), path(bed_file)
        path(bam_bai_files, stageAs: 'bam_files/')
        val(n_states)
        val(n_bins)
        val(pseudocount)

    output:
        path("${meta.model}/${meta.model}.RData"), emit: rdata
        path("${meta.model}"), emit: directory

    script:
        """
        mkdir -p ${meta.model}

        learnModel.R -r ${bed_file} -m bam_files -n ${n_states} -b ${n_bins} -p ${pseudocount} -f ${meta.model} -o ${meta.model} -t ${task.cpus}
        """
}
