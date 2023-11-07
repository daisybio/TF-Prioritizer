process LEARN_MODEL {
    tag "$meta.id"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        tuple val(meta), path(bed_file)
        tuple val(meta2), path(bam_files)
        tuple val(meta3), path(bai_files)
        val(n_states)
        val(n_bins)
        val(pseudocount)

    output:
        path("${meta.model}/${meta.model}.RData"), emit: rdata
        path("${meta.model}"), emit: directory

    script:
        """
        mkdir -p bam_files
        for file in ${bam_files}; do
            ln -fs ../\$file bam_files/\$file
        done

        for file in ${bai_files}; do
            ln -fs ../\$file bam_files/\$file
        done

        mkdir -p ${meta.model}

        learnModel.R -r ${bed_file} -m bam_files -n ${n_states} -b ${n_bins} -p ${pseudocount} -f ${meta.model} -o ${meta.model}
        """
}
