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
    create_map.py --input names.txt --taxonomy $taxonomy --output map.tsv
    """
}

process ENSG_MAPPING {
    conda "pandas"
    container "tfprio-python"

    input:
    path(count_file)
    path(map_file)

    output:
    path("ensg_counts.tsv")

    script:
    """
    convert.py --input $count_file --map $map_file --output ensg_counts.tsv
    """
}