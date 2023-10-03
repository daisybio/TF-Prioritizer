include { CREATE_ANNDATA } from "../../modules/local/counts/create_anndata"


workflow COUNTS {
    take:
        ch_counts
        ch_design

    main:
        CREATE_ANNDATA(
            ch_counts.combine(ch_design).map{counts, design -> [[id: "counts"], counts, design]},
            "gene_id"
        )
}