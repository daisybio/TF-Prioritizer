process LEARN_MODEL {

	//TODO: Check if java container needed

    input:
    path binarized_bams
    val states
    val organism

    output:
    path "ChromHMM_output/emissions_${states}.txt", emit: emissions
	path "ChromHMM_output/*_${states}_dense.bed", emit: beds

    script:
    """
    //TODO: Check for parameters for the number of states (default 10) and organism

	java -jar ${projectDir}/assets/ChromHMM.jar LearnModel \
        -p $task.cpus \
		$binarized_bams \
		ChromHMM_output \
		$states \
		$organism
    """

    stub:
	"""
	mkdir ChromHMM_output
	touch ChromHMM_output/emissions_${states}.txt
	touch ChromHMM_output/L1_${states}_dense.bed
	touch ChromHMM_output/L10_${states}_dense.bed
	touch ChromHMM_output/P6_${states}_dense.bed
	touch ChromHMM_output/P13_${states}_dense.bed
	"""
}

