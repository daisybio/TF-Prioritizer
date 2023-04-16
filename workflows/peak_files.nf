if (params.chipseq_samplesheet) { ch_chipseq_samplesheet = Channel.value(file(params.chipseq_samplesheet)) } else { exit 1, 'ChIPseq samplesheet not specified!' }

include { CHIPSEQ } from '../subworkflows/local/chipseq'
include { CREATE_FOOTPRINTS } from '../modules/local/create_footprints'
include { GUNZIP as GUNZIP_BLACKLIST } from '../subworkflows/nf-core/chipseq/modules/nf-core/modules/gunzip/main'
include { BLACKLIST } from '../modules/local/blacklist'
include { TEPIC } from '../modules/local/tepic'
include { SEQUENCE_TO_BED } from '../modules/local/sequence_to_bed'
include { MERGE_BINDING_BED } from '../modules/local/merge_binding_bed'
include { MEAN_AFFINITIES } from '../modules/local/affinities'

workflow PEAK_FILES {
    ch_versions = Channel.empty()

    if (params.chipseq_peaks) {
        ch_peaks = Channel.fromPath(params.chipseq_peaks + '/*')
    } else {
        CHIPSEQ ()
        ch_peaks = CHIPSEQ.out.peaks
        ch_versions = ch_versions.mix(CHIPSEQ.out.versions)
    }

    // TODO: Mix samples

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

    // TODO: Mix mutually exclusive

    TEPIC (
        ch_peaks, 
        params.tepic_pwm, 
        params.tepic_gtf,
        params.fasta, 
        params.tepic_windowSize,
        params.tepic_loopWindows,
        params.tepic_exponentialDecay,
        params.tepic_doNotNormalizePeakLength,
        params.tepic_maxMinutesPerChromosome,
        params.tepic_originalDecay,
        params.tepic_pValue
        )

    ch_sequences = TEPIC.out.sequences
    ch_affinities = TEPIC.out.affinities.groupTuple(by: [0, 1])
    ch_filtered_regions = TEPIC.out.filtered_regions

    ch_binding_bed_grouped = SEQUENCE_TO_BED (ch_sequences, Channel.value(params.tepic_affinityCutOff))
        .transpose(by: 2)
        .map { [it[0], it[1], it[2].name.replaceAll(/_sorted\.bed$/, ''), it[2]] }
        .groupTuple(by: [0, 1, 2])

    MERGE_BINDING_BED (ch_binding_bed_grouped)

    MEAN_AFFINITIES (ch_affinities)

    emit:
    peaks = ch_peaks
    versions = ch_versions
}