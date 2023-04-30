process CHROMOSOME_LENGTHS {
    container 'tfprio-python'

    input:
        val(species)

    output:
        path('chromosome_lengths.tsv')

    script:
        """
        get_chromosome_lengths.py -s $species -o chromosome_lengths.tsv
        """
}