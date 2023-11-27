process CONSTRUCT_MODEL {
    tag "construct"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        tuple path(background_model, stageAs: 'background.model'), path(background_counts, stageAs: 'background.counts'), path(background_bed, stageAs: 'background.bed')
        tuple path(enhancer_model, stageAs: 'enhancers.model'), path(enhancers_counts, stageAs: 'enhancers.counts'), path(enhancers_bed, stageAs: 'enhancers.bed')
        tuple path(promoter_model, stageAs: 'promoters.model'), path(promoters_counts, stageAs: 'promoters.counts'), path(promoters_bed, stageAs: 'promoters.bed')

    output:
        path('model'), emit: model_dir
        path('model/model.txt'), emit: model

    script:
        """
        mkdir -p model

        ehmm constructModel \\
            --model.e ${enhancer_model} \\
            --counts.e ${enhancers_counts} \\
            --regions.e ${enhancers_bed} \\
            --model.p ${promoter_model} \\
            --counts.p ${promoters_counts} \\
            --regions.p ${promoters_bed} \\
            --model.bg ${background_model} \\
            --nthreads ${task.cpus} \\
            --outdir model
        """
}