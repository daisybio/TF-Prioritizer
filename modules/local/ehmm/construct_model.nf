process CONSTRUCT_MODEL {
    tag "construct"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        path(background_model, stageAs: 'background.RData')
        path(enhancer_model, stageAs: 'enhancers.RData')
        path(promoter_model, stageAs: 'promoters.RData')

    output:
        path('model'), emit: model_dir
        path('model/model.RData'), emit: model

    script:
        """
        mkdir -p model

        ehmm constructModel \\
            --enhancers ${enhancer_model} \\
            --promoters ${promoter_model} \\
            --background ${background_model} \\
            --nthreads ${task.cpus} \\
            --outdir model
        """
}