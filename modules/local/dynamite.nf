process INTEGRATE_DATA {
    container 'python'

    input:
        tuple val(group1), val(group2), val(hm), path(affinity_ratios), path(diff_expression)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_integrated.tsv")

    script:
        """
        python3 $projectDir/assets/tepic/DYNAMITE/Scripts/integrateData.py \\
            $affinity_ratios \\
            $diff_expression \\
            "${group1}:${group2}_${hm}_integrated.tsv" \\
            --geneIDs 0 --expressionC 1
        """
}

process PREPARE_FOR_CLASSIFICATION {
    container 'rocker/r-ubuntu'

    input:
        tuple val(group1), val(group2), val(hm), path(integrated_data)
    
    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_classification.tsv")
    
    script:
        """
        Rscript $projectDir/assets/tepic/DYNAMITE/Scripts/prepareForClassification.R \\
            $integrated_data \\
            "${group1}:${group2}_${hm}_classification.tsv"
        """
}

process DYNAMITE {
    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    container "tepic"

    input:
        tuple val(group1), val(group2), val(hm), path(classification_data)
        val(ofolds)
        val(ifolds)
        val(alpha)
        val(randomize)
    
    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_dynamite/Regression_Coefficients_Entire_Data_Set_${group1}:${group2}_${hm}_dynamite.txt")
    
    script:
        """
        mkdir -p input
        ln -nsf ../$classification_data input/${group1}:${group2}_${hm}_dynamite.tsv

        Rscript $projectDir/assets/tepic/DYNAMITE/Scripts/DYNAMITE.R \\
            --dataDir=input \\
            --outDir="${group1}:${group2}_${hm}_dynamite" \\
            --out_var="Expression" \\
            --Ofolds=$ofolds \\
            --Ifolds=$ifolds \\
            --alpha=$alpha \\
            --performance=TRUE \\
            --randomize=$randomize
        """
}

process DYNAMITE_FILTER {
    container 'ubuntu:22.04'

    input:
        tuple val(group1), val(group2), val(hm), path(regression_coefficients)
        val(min_regression)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_dynamite_filtered.tsv")

    script:
        """
        awk -v min_regression=$min_regression 'BEGIN{OFS="\t"} NR==1{\$2="$group1:$group2"} NR==1 || \$2 >= min_regression' $regression_coefficients > ${group1}:${group2}_${hm}_dynamite_filtered.tsv
        """
}