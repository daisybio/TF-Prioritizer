process BIOPHYSICAL_MODELS {
    conda "python"
    container "tfprio-python"

    input:
        path(tfs)
        path(pwms)

    output:
        path("*.pwm"), emit: pwms
        path("*.png"), emit: plots

    script:
        """
        biophysical_models.py \
            --tfs ${tfs} \
            --pwms ${pwms}
        """
}