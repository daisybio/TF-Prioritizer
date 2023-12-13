process LEARN_MODEL {

	label "process_high"
	container "docker://openjdk:17.0.1-jdk"

    input:
    path binarized_bams
    val states

    output:
    path "ChromHMM_output/emissions_${states}.txt", emit: emissions
	path "ChromHMM_output/*_${states}_dense.bed", emit: beds

    script:
    """
	// Organism (PLACEHOLDER) only needed for downstream analysis of ChromHMM and therefore not supplied

	java -jar ${projectDir}/assets/ChromHMM.jar LearnModel \
        -p $task.cpus \
		$binarized_bams \
		ChromHMM_output \
		$states \
		PLACEHOLDER
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

