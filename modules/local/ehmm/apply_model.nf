process APPLY_MODEL {
  tag "$meta.id"
  label "process_medium"
  container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

  input:
    tuple val(meta), path(bed), val(marks), path(bam_files, stageAs: "bams/"), path(bai_files, stageAs: "bams/")
    path(model)
    path(refCounts)
    tuple val(meta2), path(index)
    val(binSize)

  output:
    tuple val(meta), path("enhancerRegions.bed"), emit: enhancers
    tuple val(meta), path("promoterRegions.bed"), emit: promoters

  script:
    marks_string = ""
    for (int i = 0; i < marks.size(); i++) {
        marks_string += " --mark ${marks[i]}:${bam_files[i]} "
    }
    """
    ehmm applyModel --regions ${bed} \\
                  --genomeSize ${index} \\
                  --model ${model} \\
                  --refCounts ${refCounts} \\
                  ${marks_string} \\
                  --nthreads ${task.cpus} \\
                  --binsize ${binSize}
    """
}
