if (params.chipseq_samplesheet) { ch_chipseq_samplesheet = Channel.value(file(params.chipseq_samplesheet)) } else { exit 1, 'ChIPseq samplesheet not specified!' }

include { CHIPSEQ } from '../subworkflows/local/chipseq'
include { CREATE_FOOTPRINTS } from '../modules/local/create_footprints'
include { GUNZIP as GUNZIP_BLACKLIST } from '../subworkflows/nf-core/chipseq/modules/nf-core/modules/gunzip/main'
include { BLACKLIST } from '../modules/local/blacklist'

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

    emit:
    peaks = ch_peaks
    versions = ch_versions
}