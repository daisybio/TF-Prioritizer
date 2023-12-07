process MAKE_CELLMARKFILETABLE {

    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    path bamDirectory

    output:
    path "cellmarkfiletable.txt"

    script:
    """
    make_cellmarkfiletable.py --input_dir $bamDirectory --output cellmarkfiletable.txt
    """

    stub:
	"""
	touch cellmarkfiletable.txt
	"""
}
