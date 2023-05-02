process EPD_NEW {
    container 'tfprio-python'

    input:
        val(genome)

    output:
        tuple val(genome), path('epd_new.bed')

    script:
        """
        fetch_epd_new.py -g ${genome} -o epd_new.bed
        sed -i 's/ /\t/g' epd_new.bed
        """
}