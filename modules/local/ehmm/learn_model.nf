process LEARN_MODEL {
    tag "$meta.id"
    label "process_medium"
    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-ehmm'

    input:
        tuple val(meta), path(bed_file)
        tuple val(marks), path(bam_files, stageAs: "bams/"), path(bai_files, stageAs: "bams/")
        val(n_states)
        val(n_bins)
        val(pseudocount)

    output:
        path("${meta.model}/model.RData"), emit: model

    script:
        marks_string = ""
        for (int i = 0; i < marks.size(); i++) {
            marks_string += " --mark ${marks[i]}:${bam_files[i]} "
        }

        """
        mkdir -p ${meta.model}

        ehmm learnModel ${marks_string} --regions ${bed_file} --nstates ${n_states} --pseudoCount ${pseudocount} --outdir ${meta.model} --nthreads ${task.cpus}
        """
}
