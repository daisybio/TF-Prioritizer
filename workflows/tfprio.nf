#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.transcript_fasta = WorkflowMain.getGenomeAttribute(params, 'transcript_fasta')
params.additional_fasta = WorkflowMain.getGenomeAttribute(params, 'additional_fasta')
params.gtf              = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff              = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed         = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.bbsplit_index    = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.star_index       = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index     = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index       = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index     = WorkflowMain.getGenomeAttribute(params, 'salmon')
params.bwa_index        = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bowtie2_index    = WorkflowMain.getGenomeAttribute(params, 'bowtie2')
params.chromap_index    = WorkflowMain.getGenomeAttribute(params, 'chromap')
params.blacklist        = WorkflowMain.getGenomeAttribute(params, 'blacklist')
params.macs_gsize       = WorkflowMain.getMacsGsize(params)
params.enable_conda     = false


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.hisat2_build_memory = params.rnaseq_hisat2_build_memory
params.fragment_size = params.chipseq_fragment_size
params.narrow_peak = params.chipseq_narrow_peak
params.min_reps_consensus = params.chipseq_min_reps_consensus
params.clip_r1 = params.chipseq_clip_r1
params.clip_r2 = params.chipseq_clip_r2
params.three_prime_clip_r1 = params.chipseq_three_prime_clip_r1
params.three_prime_clip_r2 = params.chipseq_three_prime_clip_r2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ } from '../subworkflows/local/rnaseq'
include { PEAK_FILES } from '../subworkflows/local/peak_files'
include { COLLECT_EXPRESSION } from '../modules/local/collect'
include { REPORT } from '../subworkflows/local/report'
include { DYNAMITE } from '../subworkflows/local/dynamite'
include { RANKING } from '../subworkflows/local/ranking'
include { PREPARE_TFGROUPS } from '../subworkflows/local/report'
include { CHIP_ATLAS } from '../subworkflows/local/chipatlas'
include { CHROMOSOME_LENGTHS } from '../modules/local/chromosome_lengths'
include { EHMM } from '../subworkflows/local/ehmm'

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow TFPRIO {
    ch_versions = Channel.empty()

    RNASEQ ()
    PEAK_FILES ()

    ch_versions = ch_versions.mix(RNASEQ.out.versions, PEAK_FILES.out.versions)
    ch_affinity_ratios = PEAK_FILES.out.affinity_ratios
        .map { [it[1], it[2], it[0], it[3]] } // group1, group2, hm, affinityRatios
    ch_affinity_sums = PEAK_FILES.out.affinity_sums
        .map { [it[1], it[2], it[0], it[3]] } // group1, group2, hm, affinitySums
    ch_diff_expression = RNASEQ.out.deseq2 // group1, group2, deseq2-file
    ch_counts_per_sample = RNASEQ.out.count_per_sample // group1, group2, counts-file
    ch_tpm = RNASEQ.out.tpm // group1, group2, tpm-file
    ch_counts = RNASEQ.out.count // group1, group2, counts-file
    ch_rnaseq_samplesheet = RNASEQ.out.samplesheet // samplesheet-file
    ch_peaks = PEAK_FILES.out.peaks

    ch_map = RNASEQ.out.ensg_map.splitCsv(header: false, sep: '\t')
        .map { [it[0].toUpperCase().replaceAll(/-/, '.'), it[1]]} // SYMBOL, ENSG

    DYNAMITE (
        ch_affinity_ratios,
        ch_diff_expression,
        Channel.value(params.dynamite_ofolds),
        Channel.value(params.dynamite_ifolds),
        Channel.value(params.dynamite_alpha),
        Channel.value(params.dynamite_randomize),
        Channel.value(params.dynamite_min_regression)
    )

    RANKING (
        DYNAMITE.out,
        ch_diff_expression,
        ch_affinity_sums,
        ch_map
    )

    PREPARE_TFGROUPS (
        RANKING.out.groups,
        Channel.value(params.top_target_genes),
        ch_affinity_sums,
        ch_counts_per_sample,
        ch_rnaseq_samplesheet,
        Channel.value(params.tepic_pwm),
        RNASEQ.out.ensg_map
    )

    COLLECT_EXPRESSION (
        ch_tpm,
        RNASEQ.out.deseq2_grouped,
        ch_counts,
        RANKING.out.ensgs
    )

    if (params.chipatlas_tissue_types != null)
    { 
        chr_lengths = CHROMOSOME_LENGTHS(Channel.value(params.biomart_species))

        ch_chipatlas = CHIP_ATLAS(
            Channel.value(params.chipatlas_tissue_types),
            RANKING.out.tfs.splitCsv(header: false).map{ it[0].toLowerCase() }.collect(),
            Channel.value(params.chipatlas_genome),
            Channel.value(params.chipatlas_threshold),
            chr_lengths
        )

        EHMM(
            ch_chipatlas,
            chr_lengths,
            Channel.value(params.chipatlas_genome),
            Channel.value(params.chipatlas_tissue_types),
            ch_peaks
        )
    }

    /*REPORT (
        COLLECT_EXPRESSION.out,
        RANKING.out.ranks,
        PREPARE_TFGROUPS.out.flatten().collect(),
        RANKING.out.group_tfs_map,
        RANKING.out.tf_ensg_map
    )*/
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    TFPRIO ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/