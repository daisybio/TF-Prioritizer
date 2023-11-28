include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"
include { APPLY_MODEL } from "../../modules/local/ehmm/apply_model"
include { CAT_CAT as CAT_BACKGROUND } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as MERGE_ASSAYS } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as MERGE_PROMOTERS_ENHANCERS } from "../../modules/nf-core/cat/cat/main"

// Groovy function mapping a directory to a list of files in it with a matching extension
def directory_to_files(dir, ext) {
    def files = []
    dir.eachFileRecurse () { file ->
        if (file.name.endsWith(ext)) {
            files.add(file)
        }
    }
    return files
}

workflow EHMM {
    take:
        background
        enhancers_bed
        promoters_bed
        peaks
        bams

        gtf
        index

        genomic_region_size
        train_split
        random_seed
        top_quantile
        n_samples

        n_states
        n_bins
        pseudocount

    main:
        CAT_BACKGROUND(
            background.map{ it[1] }.collect().map { [[id: 'background'], it] }
        )

        SPLIT_DATASETS(
            CAT_BACKGROUND.out.file_out,
            enhancers_bed,
            promoters_bed,
            gtf,
            genomic_region_size,
            train_split,
            random_seed,
            top_quantile,
            n_samples
        )

        ch_bam_bai = background.map{ meta, bed, bam, bai -> [meta["id"], bam, bai]}
                               .reduce([[], [], []]){ acc, curr -> [acc[0] + [curr[0]], acc[1] + [curr[1]], acc[2] + [curr[2]]]}

        LEARN_BACKGROUND_MODEL (
            SPLIT_DATASETS.out.trainBackground.map{[[id: "train_background", model: "BackgroundModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers.map{[[id: "train_enhancers", model: "EnhancerModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters.map{[[id: "train_promoters", model: "PromoterModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        CONSTRUCT_MODEL (
            LEARN_BACKGROUND_MODEL.out.model,
            LEARN_ENHANCER_MODEL.out.model,
            LEARN_PROMOTER_MODEL.out.model
        )

        MERGE_ASSAYS (
            peaks.map{meta, file_ -> [[id: meta["state"], state: meta["state"]], file_]}.groupTuple()
        )

        ch_bams = bams.map{ meta, directory -> [meta["state"], meta["assay"], directory]}
                      .map{ state, assay, directory -> [state, assay, directory_to_files(directory, ".bam"), directory_to_files(directory, ".bai")]}
                      .transpose(by: [2, 3])
                      .groupTuple()
                      .map{ state, assays, bams, bais -> [[id: state, state: state], assays, bams, bais]}

        ch_combined = MERGE_ASSAYS.out.file_out .map{meta, bed -> [meta["state"], bed]}
                                                .join(ch_bams.map{meta, assays, bams, bais -> [meta["state"], assays, bams, bais]})
                                                .map{state, bed, assays, bams, bais -> [[id: state], bed, assays, bams, bais]}

        APPLY_MODEL (
            ch_combined,
            CONSTRUCT_MODEL.out.model,
            index,
            n_bins
        )

        ch_enhancers = APPLY_MODEL.out.enhancers
                        .map{meta, regions -> [[id: meta["state"] + "_enhancers", 
                                                state: meta["state"], 
                                                antibody: "enhancers"], regions]}
                        .filter{meta, regions -> regions.size() > 0}
        ch_promoters = APPLY_MODEL.out.promoters
                        .map{meta, regions -> [[id: meta["state"] + "_promoters", 
                                                state: meta["state"], 
                                                antibody: "promoters"], regions]}
                        .filter{meta, regions -> regions.size() > 0}

        /*
        MERGE_PROMOTERS_ENHANCERS (
            ch_enhancers.map{meta, regions -> [[id: meta["state"]], regions]}.mix(
                ch_promoters.map{meta, regions -> [[id: meta["state"]], regions]}
            ).groupTuple()
        )

        ch_all = MERGE_PROMOTERS_ENHANCERS.out.file_out.map{meta, regions -> [[id: meta["id"] + "_all", 
                                                                                state: meta["id"], 
                                                                                antibody: "regulatoryElements"], regions]}
        */

        ch_all = ch_enhancers.mix(ch_promoters)
    
    emit:
        enhancers = ch_enhancers
        promoters = ch_promoters
        all = ch_all
}
