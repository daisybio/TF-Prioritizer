process REFORMAT_BAM {

    label "process_medium"
    container "registry.hub.docker.com/staphb/samtools"

    input:
    path bamFileIn

    output:
    path "${bamFileIn.baseName}_reformatted.bam"

    script:
    """
    reformat_bam.sh $bamFileIn ${bamFileIn.baseName}_reformatted.bam
    """

	stub:
	"""
	touch ${bamFileIn.baseName}_reformatted.bam
	"""
}

