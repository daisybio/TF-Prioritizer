process LEARN_MODEL {

	//TODO: Check if java container needed

    input:
    path binarized_bams

    output:
    path "ChromHMM_output"

    script:
    """
    //TODO: Check for parameters for the number of states (default 10) and organism

	java -jar ${projectDir}/assets/ChromHMM.jar LearnModel \
		$binarized_bams \
		ChromHMM_output \
		10 \
		mm10
    """
}

