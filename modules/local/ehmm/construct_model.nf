process CONSTRUCT_MODEL {
    container 'tfprio-r'
    label 'process_medium'

    input:
        path(background_model)
        path(enhancer_model)
        path(promoter_model)

    output:
        path('model'), emit: model_dir
        path('model/model.txt'), emit: model

    script:
        """
        mkdir model

        constructModel.R \\
            -b ${background_model} \\
            -e ${enhancer_model} \\
            -p ${promoter_model} \\
            -t ${task.cpus} \\
            -o model
        """
}