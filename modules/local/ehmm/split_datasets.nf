process SPLIT_DATASETS {
    tag "split_datasets"
    label 'process_high_memory'

    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        tuple val(meta), path(background)
        tuple val(meta2), path(enhancers)
        tuple val(meta3), path(promoters)
        path(gtf)
        val(genomic_region_size)
        val(train_split)
        val(random_seed)
        val(top_quantile)
        val(sample_size)

    output:
        path('testBackground.bed'), emit: testBackground
        path('testEnhancers.bed'), emit: testEnhancers
        path('testPromoters.bed'), emit: testPromoters
        path('trainBackground.bed'), emit: trainBackground
        path('trainEnhancers.bed'), emit: trainEnhancers
        path('trainPromoters.bed'), emit: trainPromoters

    script:
        """
        split_datasets.R \\
            -b ${background} \\
            -e ${enhancers} \\
            -p ${promoters} \\
            -g ${gtf} \\
            -o . \\
            -f ${sample_size} \\
            -s ${genomic_region_size} \\
            -r ${random_seed} \\
            -q ${top_quantile} \\
            -t ${train_split}
        """
}
