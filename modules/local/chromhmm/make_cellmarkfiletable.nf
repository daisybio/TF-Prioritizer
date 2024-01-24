process MAKE_CELLMARKFILETABLE {

    label "process_single"
    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    path bam_design

    output:
    path "cellmarkfiletable.txt"

    script:
    """
    make_cellmarkfiletable.py --input $bam_design --output cellmarkfiletable.txt
    """

    stub:
	"""
	touch cellmarkfiletable.txt
	"""
}
