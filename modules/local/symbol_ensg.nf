process SYMBOL_ENSG {
    container "tfprio-python"

    input:
        path(symbol_ensg)

    output:
        path("symbol_ensg.json")

    script:
        """
        df_to_json.py -i $symbol_ensg -o symbol_ensg.json
        """
}