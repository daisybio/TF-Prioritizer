include { CREATE_ANNDATA } from "../../modules/local/counts/create_anndata"
include { FILTER_ANNDATA } from "../../modules/local/counts/filter_anndata"
include { FETCH_SYMBOL_MAP } from "../../modules/local/counts/fetch_symbol_map"
include { USE_SYMBOL_MAP } from "../../modules/local/counts/use_symbol_map"

workflow COUNTS {
    take:
        ch_counts
        ch_design

    main:
        CREATE_ANNDATA(
            ch_counts.combine(ch_design).map{counts, design -> [[id: "counts"], counts, design]},
            "gene_id"
        )

        FILTER_ANNDATA(
            CREATE_ANNDATA.out.anndata,
            params.min_count,
            params.min_tpm
        )

        FETCH_SYMBOL_MAP(
            CREATE_ANNDATA.out.genes,
            params.tax_id
        )

        USE_SYMBOL_MAP(
            FILTER_ANNDATA.out,
            FETCH_SYMBOL_MAP.out
        )
}