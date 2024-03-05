process RUN_FIMO {

    container "registry.hub.docker.com/leonhafner/memesuite"

    input:
        path motif_file
        path sequence_file
    
    output:
        path "fimo_out"
    
    script:
    """
    fimo ${motif_file} ${sequence_file}
    """

    stub:
    """
    mkdir fimo_out
    touch fimo_out/best_site.narrowPeak
    touch fimo_out/cisml.xml
    touch fimo_out/fimo.gff
    touch fimo_out/fimo.html
    touch fimo_out/fimo.tsv
    touch fimo_out/fimo.xml
    """

    
}