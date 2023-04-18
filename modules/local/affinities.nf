process MEAN_AFFINITIES {
    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
        tuple val(hm), val(group), path(affinities)

    output:
        tuple val(hm), val(group), path("${hm}_${group}_meanAffinities.txt")

    script:
    if (affinities.getClass().isArray() && affinities.size() > 1) {
        """
        mean_affinities.R --input ${affinities.join(' ')} --output ${hm}_${group}_meanAffinities.txt
        """
    }
    else {
        """
        ln -s ${affinities[0]} ${hm}_${group}_meanAffinities.txt
        """
    }
}

process AFFINITY_RATIOS {
    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
        tuple val(hm), val(group1), path(affinities1), val(group2), path(affinities2)

    output:
        tuple val(hm), val(group1), val(group2), path("${hm}_${group1}:${group2}_affinityRatios.txt")

    script:
    """
    affinity_ratios.R --input1 ${affinities1} --input2 ${affinities2} --output ${hm}_${group1}:${group2}_affinityRatios.txt
    """
}