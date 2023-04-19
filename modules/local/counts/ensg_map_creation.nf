process ENSG_MAP_CREATION {
    conda "pandas mygene"
    container "tfprio-python"

    input:
    path(count_file)
    val(taxonomy)

    output:
    path("map.tsv")

    script:
    """
    cut -f 1 ${count_file} | sed '1d' | sort -u > names.txt
    convert.py --input names.txt --taxonomy $taxonomy --output map.tsv
    """
}