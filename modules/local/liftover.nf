process LIFTOVER {
    tag "$meta.id"
    label 'process_single'

    container 'registry.hub.docker.com/bigdatainbiomedicine/easyliftover:latest'

    input:
        tuple val(meta), path(input)
        val current_genome
        val target_genome

    output:
        tuple val(meta), path("lifted_${input.name}")

    script:
        """
        #!/usr/bin/env python

        from easyliftover import liftover_path

        result = liftover_path("$current_genome", "$target_genome", "$input")

        with open("lifted_${input.name}", "w+") as f:
            f.write(result)
        """
}
