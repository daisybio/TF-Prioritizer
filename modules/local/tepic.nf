process TEPIC {
    conda "bioconda::bedtools conda-forge::gxx conda-forge::r-base conda-forge::r-gplots conda-forge::r-ggplot2 conda-forge::r-glmnet conda-forge::r-doMC conda-forge::r-reshape2 conda-forge::r-gridExtra pandas"
    container "tepic"
    cpus 4

    input:
    tuple val(sample), path(experimental_bed), val(hm), val(group)
    path(pwm)
    path(gtf)
    path(reference_fasta)
    val(windowSize)
    val(loopWindows)
    val(exponentialDecay)
    val(doNotNormalizePeakLength)
    val(maxMinutesPerChromosome)
    val(originalDecay)
    val(pValue)

    output:
    path("sequences.tsv"), emit: sequences
    path("regions_to_target_genes.tsv"), emit: regions_to_target_genes
    path("*_Affinity.txt"), emit: affinity
    path("*_canidate_binding_regions_Filtered_Regions.bed"), emit: filtered_regions
    path("*_Decay_Peak_Features_Affinity_Gene_View.txt"), emit: gene_view

    // TODO: Implement missing config parameters

    script:
    """
    g++ $projectDir/assets/tepic/code/TRAPmulti.cpp -O3 -fopenmp -o $projectDir/assets/tepic/code/TRAPmulti

    bash $projectDir/assets/tepic/code/TEPIC.sh \
        -b $experimental_bed \
        -o $sample \
        -S sequences.tsv \
        -g $reference_fasta \
        -a $gtf \
        -p $pwm \
        -w $windowSize \
        -s $loopWindows \
        -e $exponentialDecay \
        -l $doNotNormalizePeakLength \
        -i $maxMinutesPerChromosome \
        -x $originalDecay \
        -v $pValue \
        -j \
        -c $task.cpus
    """
}