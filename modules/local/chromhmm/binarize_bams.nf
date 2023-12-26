process BINARIZE_BAMS {

	label "process_high"
	container "registry.hub.docker.com/leonhafner/openjdk:17"

    input:
    path bams, stageAs: "reformatted_bams/*"
	path bais, stageAs: "reformatted_bams/*"
	path cellmarkfiletable
	path chromsizes

    output:
    path "binarized_bams"

    script:
    """
	ChromHMM.jar BinarizeBam \
		$chromsizes \
		reformatted_bams \
		$cellmarkfiletable \
		binarized_bams
	"""

	stub:
	"""
	mkdir binarized_bams
	touch binarized_bams/chr1.txt
	touch binarized_bams/chr2.txt
	"""
}
