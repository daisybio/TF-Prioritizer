include { MAKE_CELLMARKFILETABLE } from "../../modules/local/chromhmm/make_cellmarkfiletable"
include { REFORMAT_BAM } from "../../modules/local/chromhmm/reformat_bam"
include { SAMTOOLS_INDEX as INDEX_BAM } from "../../modules/nf-core/samtools/index/main"
include { BINARIZE_BAMS } from "../../modules/local/chromhmm/binarize_bams"
include { LEARN_MODEL } from "../../modules/local/chromhmm/learn_model"
include { GET_OUTPUT } from "../../modules/local/chromhmm/get_output"



workflow CHROMHMM {
	take:
		raw_bams
		chromhmm_states // default to 10
		chromsizes


	main:
		ch_bams_raw = Channel.fromPath("${raw_bams}/*/*.bam")

    	MAKE_CELLMARKFILETABLE("${raw_bams}")

    	REFORMAT_BAM(ch_bams_raw)

    	INDEX_BAM(REFORMAT_BAM.out)

    	BINARIZE_BAMS(INDEX_BAM.out.bai.collect(), MAKE_CELLMARKFILETABLE.out, chromsizes)

    	LEARN_MODEL(BINARIZE_BAMS.out, chromhmm_states)

		ch_emission_bed = LEARN_MODEL.out.emissions.combine(LEARN_MODEL.out.beds.flatten())

		GET_OUTPUT(ch_emission_bed)
}
