process FILTER_MOTIFS {
    container "registry.hub.docker.com/bigdatainbiomedicine/inspect-python"

    input:
        path(tfs_jaspar_ids)
        path(jaspar_motifs)
    
    output:
        path "motifs.meme"
    
    script:
    """
        #!/usr/bin/env python3

        from os.path import exists

        tfs_jaspar_ids = "${tfs_jaspar_ids}"
        jaspar_motifs = "${jaspar_motifs}"

        # Read differentially expressed (DE) transcription factors (TF)
        with open(tfs_jaspar_ids, "r") as f:
            tfs_jaspar_ids = f.read().split('\\n')

        meme_motifs_file = []

        # Iterate over TFs and store meme files for DE TFs
        for jaspar_id in tfs_jaspar_ids:
            if exists(f"jaspar_motifs/{jaspar_id}.meme"):
                with open(f"jaspar_motifs/{jaspar_id}.meme", "r") as f:
                    # Skip first 9 lines, as they are the same for every file
                    lines = "\\n".join(f.read().split("\\n")[9:])
                    meme_motifs_file.append(lines)

        # Create header for final file
        prefix = "MEME version 4\\n\\nALPHABET= ACGT\\n\\nstrands: + -\\n\\nBackground letter frequencies\\nA 0.25 C 0.25 G 0.25 T 0.25\\n\\n\\n"

        # Join header with concatenated meme files
        with open("motifs.meme", "w") as f:
            f.write(prefix + "\\n".join(meme_motifs_file))
    """

    stub:
    """
    touch motifs.meme
    """
}