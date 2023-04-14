process TEPIC {
    conda "conda-forge::r-base"
    container "r-base"

    input:
    tuple val(sample), path(experimental_bed), val(hm), val(group)
    path(pwm)
    path(reference_fasta)
    val(windowSize)
    val(loopWindows)
    val(exponentialDecay)
    val(doNotNormalizePeakLength)
    val(maxMinutesPerChromosome)
    val(originalDecay)
    val(pValue)

    script:
    """
    tepic/compile_TRAP_install_R_packages.sh

    tepic/TEPIC.sh \
        -b $experimental_bed \
        -a $pwm \
        -g $reference_fasta \
        -w $windowSize \
        -s $loopWindows \
        -e $exponentialDecay \
        -l $doNotNormalizePeakLength \
        -i $maxMinutesPerChromosome \
        -x $originalDecay \
        -v $pValue \
        -c $process.cpus
    """
}