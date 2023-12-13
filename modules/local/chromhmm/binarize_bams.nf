process BINARIZE_BAMS {

	label "process_high"
	container "docker://openjdk:17.0.1-jdk"

    input:
    path cellmarkfiletable
    path bams, stageAs: "reformatted_bams/*"
	path chromsizes

    output:
    path "binarized_bams"

    script:
    """
	java -jar ChromHMM.jar BinarizeBam \
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
