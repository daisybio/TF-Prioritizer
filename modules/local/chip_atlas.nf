process CHIP_ATLAS {
    container 'tfprio-python'

    input:
        val(tissue_types)
        val(tfs)
        val(genome)

    output:
        path('*.bed')

    script:
        """
        if [ ! -f chip_atlas_file_list.csv ]; then
            curl -o overview.zip http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv
            unzip overview.zip
        fi

        chipatlas.py -l chip_atlas_file_list.csv -t ${tissue_types.join(' ')} -f ${tfs.join(' ')} -g ${genome}
        """
}