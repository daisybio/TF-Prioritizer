process APPLY_MODEL {
  tag "$meta.id"
  label "process_medium"
  container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

  input:
    tuple val(meta), path(bed)
    path(bam_bai_files, stageAs: 'bam_files/')
    path(model)
    tuple val(meta2), path(index)
    val(binSize)
    val(pseudocount)
  
  output:
  
  script:
  """
  applyModel.R -r ${bed} \\
                -g ${index} \\
                -m ${model} \\
                -b bam_files \\
                -t ${task.cpus} \\
                -s ${binSize} \\
                -p ${pseudocount}
  """
}