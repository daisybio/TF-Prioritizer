process EPD_NEW {
    container 'tfprio-python'

    input:
        val(genome)

    output:
        path('epd_new.bed')

    script:
        """
        fetch_epd_new.py -g ${genome} -o epd_new.bed
        """
}