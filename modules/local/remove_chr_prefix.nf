process REMOVE_CHR_PREFIX {
    container 'ubuntu'

    input:
        tuple val(meta), path(bed_file)

    output:
        tuple val(meta), path("${bed_file.name.replaceFirst(/\.bed$/, '.chr.bed')}")

    script:
        """
        sed -e 's/^chr//' ${bed_file} > ${bed_file.name.replaceFirst(/\.bed$/, '.chr.bed')}
        """
}