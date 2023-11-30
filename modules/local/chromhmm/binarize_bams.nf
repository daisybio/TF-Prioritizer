process BINARIZE_BAMS {

	//TODO: Check if java container needed

    input:
    path cellmarkfiletable
    path bams, stageAs: "reformatted_bams/*"

    output:
    path "binarized_bams"

    script:
    """
	//TODO: Get parameters for organism
	
	java -jar ${projectDir}/assets/ChromHMM.jar BinarizeBam \
		${projectDir}/assets/CHROMSIZES/mm10.txt \
		reformatted_bams \
		$cellmarkfiletable \
		binarized_bams
    
	"""
}
