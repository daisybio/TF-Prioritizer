process SPLIT_DATASETS {
    container 'tfprio-r'

    input:
        path(background)
        path(enhancers)
        path(promoters)
        path(gtf)
        val(genomic_region_size)
        val(train_split)
        val(random_seed)
        val(top_quantile)
        val(sample_size)

    output:
        path('test.txt')
    
    script:
        """
        split_datasets.R \\
            -b ${background} \\
            -e ${enhancers} \\
            -p ${promoters} \\
            -g ${gtf} \\
            -o test \\
            -f ${sample_size} \\
            -s ${genomic_region_size} \\
            -r ${random_seed} \\
            -q ${top_quantile} \\
            -t ${train_split}
        """
}