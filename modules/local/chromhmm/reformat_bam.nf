process REFORMAT_BAM {

    label "process_medium"
    container "registry.hub.docker.com/staphb/samtools"

    input:
    tuple val(meta), path(bamFileIn)

    output:
    tuple val(meta), path("${bamFileIn.baseName}_reformatted.bam")

    script:
    """
    reformat_bam.sh $bamFileIn ${bamFileIn.baseName}_reformatted.bam
    """

	stub:
	"""
	touch ${bamFileIn.baseName}_reformatted.bam
	"""
}

