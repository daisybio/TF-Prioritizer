process CONSTRUCT_MODEL {
    tag "construct"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        path(background_model)
        path(enhancer_model)
        path(promoter_model)

    output:
        path('model'), emit: model_dir
        path('model/model.txt'), emit: model

    script:
        """
        mkdir -p model

        constructModel.R \\
            -b ${background_model} \\
            -e ${enhancer_model} \\
            -p ${promoter_model} \\
            -t ${task.cpus} \\
            -o model
        """
}