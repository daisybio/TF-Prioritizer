include { JASPAR_MAPPING } from "../../modules/local/fimo/jaspar_mapping"
include { JASPAR_DOWNLOAD } from "../../modules/local/fimo/jaspar_download"
include { FILTER_MOTIFS } from "../../modules/local/fimo/filter_motifs"
include { CAT_CAT as CONCAT_BEDS} from "../../modules/nf-core/cat/cat/main"
include { BEDTOOLS_SORT as SORT_REGIONS } from "../../modules/nf-core/bedtools/sort/main"
include { BEDTOOLS_MERGE as MERGE_REGIONS } from "../../modules/nf-core/bedtools/merge/main"
include { BEDTOOLS_GETFASTA as EXTRACT_SEQUENCE } from "../../modules/nf-core/bedtools/getfasta/main"
include { RUN_FIMO } from "../../modules/local/fimo/run_fimo"


workflow FIMO {
	take:
		tfs_sorted
        enhancer_regions

	main:

		JASPAR_MAPPING(tfs_sorted)

		JASPAR_DOWNLOAD()

		FILTER_MOTIFS(JASPAR_MAPPING.out, JASPAR_DOWNLOAD.out)
	
		ch_cat_input = enhancer_regions
			.map{
				meta, file -> file
			}
			.collect()
			.map{
				item -> [[id: "enhancer_regions"], item]
			}
		
		CONCAT_BEDS(ch_cat_input)

		SORT_REGIONS(CONCAT_BEDS.out.file_out, [])

		MERGE_REGIONS(SORT_REGIONS.out.sorted)

		ch_bed = MERGE_REGIONS.out.bed.map{meta, bed -> bed}

		EXTRACT_SEQUENCE(ch_bed, params.fasta)

		RUN_FIMO(FILTER_MOTIFS.out, EXTRACT_SEQUENCE.out.fasta)
		RUN_FIMO.out.view()
}
