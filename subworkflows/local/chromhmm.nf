include { MAKE_CELLMARKFILETABLE } from "../../modules/local/chromhmm/make_cellmarkfiletable"
include { REFORMAT_BAM } from "../../modules/local/chromhmm/reformat_bam"
include { SAMTOOLS_INDEX as INDEX_BAM } from "../../modules/nf-core/samtools/index/main"
include { BINARIZE_BAMS } from "../../modules/local/chromhmm/binarize_bams"
include { LEARN_MODEL } from "../../modules/local/chromhmm/learn_model"



workflow CHROMHMM {
	take:
		raw_bams

	main:
		ch_bams_raw = Channel.fromPath("${raw_bams}/*/*.bam")

    	MAKE_CELLMARKFILETABLE("${raw_bams}")

    	REFORMAT_BAM(ch_bams_raw)

    	INDEX_BAM(REFORMAT_BAM.out)

    	BINARIZE_BAMS(MAKE_CELLMARKFILETABLE.out, INDEX_BAM.out.bai.collect())

    	LEARN_MODEL(BINARIZE_BAMS.out)
}
