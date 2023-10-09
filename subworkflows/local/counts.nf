include { CREATE_DATAFRAME } from "../../modules/local/counts/create_dataframe"
include { FILTER_ANNDATA } from "../../modules/local/counts/filter_anndata"
include { FETCH_SYMBOL_MAP } from "../../modules/local/counts/fetch_symbol_map"
include { CREATE_ANNDATA } from "../../modules/local/counts/create_anndata"
include { USE_SYMBOL_MAP } from "../../modules/local/counts/use_symbol_map"
include { NORMALIZE } from "../../modules/local/counts/normalize"
include { DESEQ2 } from "../../modules/local/counts/deseq2"

workflow COUNTS {
    take:
        ch_counts
        ch_design

    main:
        CREATE_DATAFRAME(
            ch_counts.combine(ch_design).map{counts, design -> [[id: "counts"], counts, design]},
            "gene_id"
        )

        FETCH_SYMBOL_MAP(
            CREATE_DATAFRAME.out.genes,
            params.tax_id
        )

        USE_SYMBOL_MAP(
            CREATE_DATAFRAME.out.dataframe,
            FETCH_SYMBOL_MAP.out
        )

        CREATE_ANNDATA(
            USE_SYMBOL_MAP.out.combine(ch_design)
        )

        FILTER_ANNDATA(
            CREATE_ANNDATA.out,
            params.min_count,
            params.min_tpm
        )

        DESEQ2(
            FILTER_ANNDATA.out
        )

        ch_deseq2_stats = DESEQ2.out.transpose().map{ meta, stats_file -> 
                base = stats_file.baseName.split(":")
                state1 = base[0]
                state2 = base[1]
                [[id: state1 + ":" + state2, state1: state1, state2: state2], stats_file]
            }
    
    emit:
        deseq2 = ch_deseq2_stats
        adata = FILTER_ANNDATA.out.collect()
}