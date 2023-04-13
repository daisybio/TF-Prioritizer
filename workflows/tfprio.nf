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

if (params.rnaseq_samplesheet) { ch_rnaseq_samplesheet = file(params.rnaseq_samplesheet) } else { exit 1, 'RNAseq samplesheet not specified!' }
if (params.chipseq_samplesheet) { ch_chipseq_samplesheet = Channel.fromPath(params.chipseq_samplesheet) } else { exit 1, 'ChIPseq samplesheet not specified!' }

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
include { CHIPSEQ } from '../subworkflows/local/chipseq'
include { COUNT_PREPROCESSING } from '../modules/local/count_preprocessing'
include { CREATE_FOOTPRINTS } from '../modules/local/create_footprints'
include { GUNZIP as GUNZIP_BLACKLIST } from '../subworkflows/nf-core/chipseq/modules/nf-core/modules/gunzip/main'
include { BLACKLIST } from '../modules/local/blacklist'

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow TFPRIO {
    ch_versions = Channel.empty()

    if (params.rnaseq_counts) {
        ch_count = file(params.rnaseq_counts)
    } else {
        RNASEQ ()
        ch_count = RNASEQ.out.counts
        ch_bigwig = RNASEQ.out.bigwig
        ch_versions = ch_versions.mix(RNASEQ.out.versions)
    }

    if (params.chipseq_peaks) {
        ch_peaks = Channel.fromPath(params.chipseq_peaks + '/*')
    } else {
        CHIPSEQ ()
        ch_peaks = CHIPSEQ.out.peaks
        ch_versions = ch_versions.mix(CHIPSEQ.out.versions)
    }

    ch_peaks = CREATE_FOOTPRINTS (params.peakBindingSiteSearch, params.maxDistance, ch_peaks).footprints
        .map { [it.name, it] }
        .map { [it[0].replaceAll(/_footprints.bed$/,''), it[1]] }
        .map { [it[0].replaceAll(/_peaks$/, ''), it[1]] }

    ch_chipseq_annotations = ch_chipseq_samplesheet
        .splitCsv(header: true, sep: ',')
        // Remove control samples
        .filter { it.antibody != '' && it.group != '' }
        .map { [it.sample, it.antibody, it.group ] }
        .unique()
    
    ch_peaks = ch_peaks.combine(ch_chipseq_annotations, by: 0)

    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], params.blacklist ] ).gunzip.map{ it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.value(file(params.blacklist))
        }
    }

    ch_peaks = BLACKLIST (ch_blacklist, ch_peaks)

    // COUNT_PREPROCESSING (ch_count, ch_rnaseq_samplesheet)
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