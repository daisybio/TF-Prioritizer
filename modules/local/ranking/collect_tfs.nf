process COLLECT_TFS {
    label 'process_single'

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
        path(rankings)
    
    output:
        path("tfgroup_tfs.json"), emit: tfgroup_tfs
        path("tfs_sorted.txt"), emit: tfs
    
    script:
    // Case that rankings contains only a single element
    if (!(rankings instanceof List)) {
        rankings = [rankings]
    }
    """
    #!/usr/bin/env python3

    import json

    paths = ["${rankings.name.join("\", \"")}"]
    tfgroup_tfs = {}

    for path in paths:
        with open(path) as f:
            for line in f:
                tfgroup = line.strip().split("\\t")[0]
                tfs = tfgroup.split("..")
                tfgroup_tfs[tfgroup] = tfs

    with open("tfgroup_tfs.json", "w") as f:
        json.dump(tfgroup_tfs, f, indent=4)

    tfs_sorted = sorted(set([tf for tfs in tfgroup_tfs.values() for tf in tfs]))

    with open("tfs_sorted.txt", "w") as f:
        f.write("\\n".join(tfs_sorted))
    """
}