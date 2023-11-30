process REFORMAT_BAM {

    container "registry.hub.docker.com/staphb/samtools"

    input:
    path bamFileIn

    output:
    path "${bamFileIn.baseName}_reformatted.bam"

    script:
    """
    reformat_bam.sh $bamFileIn ${bamFileIn.baseName}_reformatted.bam
    """
}

