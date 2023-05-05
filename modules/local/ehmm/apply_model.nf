process APPLY_MODEL {
    container 'tfprio-r'
    
    input:
        tuple val(meta), path(regions)
        path(chromosome_lengths)
        path(model)
        val(bin_size)
        val(pseudocount)

    output:

    script:
        """
        applyModel.R \\
            -r ${regions} \\
            -g ${chromosome_lengths} \\
            -m ${model} \\
            
            -t ${task.cpus} \\
            -p ${pseudocount} \\
            -s ${bin_size}
        """
}