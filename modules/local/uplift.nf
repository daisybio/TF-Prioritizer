process UPLIFT {
    container 'tfprio-python'

    input:
        tuple val(meta), path(bed_file)
        val(genome)

    output:
        tuple val(meta), path("*_uplifted.bed")

    script:
        """
        uplift.py ${bed_file} ${meta.id}_uplifted.bed ${meta.genome} $genome 0 1 2
        """
}