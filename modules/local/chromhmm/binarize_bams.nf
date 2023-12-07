process BINARIZE_BAMS {

	//TODO: Check if java container needed

    input:
    path cellmarkfiletable
    path bams, stageAs: "reformatted_bams/*"
	val organism

    output:
    path "binarized_bams"

    script:
    """
	java -jar ${projectDir}/assets/ChromHMM.jar BinarizeBam \
		${projectDir}/assets/CHROMSIZES/"${organism}".txt \
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
