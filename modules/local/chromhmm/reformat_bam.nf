process REFORMAT_BAM {

    label "process_medium"
    container "registry.hub.docker.com/staphb/samtools"

    input:
    tuple val(meta), path(bamFileIn, stageAs: "input/*")
    // stage input to avoid name clashes

    output:
    tuple val(meta), path("${bamFileIn.baseName}.${bamFileIn.extension}")

    script:
    """
    reformat_bam.sh $bamFileIn ${bamFileIn.baseName}.${bamFileIn.extension}
    """

	stub:
	"""
	touch ${bamFileIn.baseName}.${bamFileIn.extension}
	"""
}

