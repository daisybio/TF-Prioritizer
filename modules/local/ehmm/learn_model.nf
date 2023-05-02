process LEARN_MODEL {
    container 'tfprio-r'

    input:
        path(bed_file)
        path(bam_files)
        path(bai_files)
        val(n_states)
        val(n_bins)
        val(pseudocount)

    output:
        path("${task.ext.model}/${task.ext.model}.RData"), emit: rdata
        path("${task.ext.model}"), emit: directory

    script:
        """
        mkdir -p bam_files
        for file in ${bam_files}; do
            ln -s ../\$file bam_files/\$file
        done

        for file in ${bai_files}; do
            ln -s ../\$file bam_files/\$file
        done

        mkdir -p ${task.ext.model}

        learnModel.R -r ${bed_file} -m bam_files -n ${n_states} -b ${n_bins} -p ${pseudocount} -f ${task.ext.model} -o ${task.ext.model}
        """
}