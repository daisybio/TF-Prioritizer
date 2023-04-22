process HEATMAPS {
    conda "conda-forge::r-base conda-forge::r-optparse conda-forge::r-pheatmap"
    container "tfprio-r"

    input:
        tuple val(group1), val(group2), val(hm), path(top_target_genes)
        path(counts)
        path(samplesheet)
        path(ensg_map)

    output:
        tuple val(group1), val(group2), val(hm), path("*.png")

    script:
        """
        heatmap.R --counts ${counts} --samplesheet ${samplesheet} --target_genes ${top_target_genes} --group1 ${group1} --group2 ${group2} --hm ${hm} --map ${ensg_map}
        """
}