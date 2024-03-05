process JASPAR_MAPPING {
    label 'process_single'

    container "registry.hub.docker.com/bigdatainbiomedicine/inspect-python"
    
    input:
        path(tfs_sorted)

    output:
        path("tfs_jaspar_ids.txt")

    script:
    """
        #!/usr/bin/env python3

        from collections import defaultdict
        from urllib.request import urlopen

        path_tfs_sorted = "${tfs_sorted}"

        # Read differential expressed TFs
        with open(path_tfs_sorted, 'r') as f:
            tfs_sorted = f.read().split("\\n")

        # Get mapping file 1
        fetch1 = urlopen("https://raw.githubusercontent.com/SchulzLab/STARE/main/PWMs/2.2/Jaspar_Hocomoco_Kellis_human_PSEMs.txt")
        mapping1 = [tuple(line[1:].split("\\t")[:2]) for line in fetch1.read().decode('utf-8').split('\\n') if line.startswith('>')]

        # Get mapping file 2
        fetch2 = urlopen("https://raw.githubusercontent.com/SchulzLab/STARE/main/PWMs/2.2/Jaspar_Hocomoco_Kellis_mouse_PSEMs.txt")
        mapping2 = [tuple(line[1:].split("\\t")[:2]) for line in fetch2.read().decode('utf-8').split('\\n') if line.startswith('>')]

        # Create mapping dict from mapping files
        symbol_to_id = defaultdict(set)
        for jaspar_id, symbol in mapping1 + mapping2:
            symbol_to_id[symbol].add(jaspar_id)

        # Cast defaultdict to dict
        symbol_to_id = dict(symbol_to_id)

        # Create file with sorted TF meme IDs
        tfs = sorted([jaspar_id for tf in tfs_sorted if tf in symbol_to_id for jaspar_id in symbol_to_id[tf]])

        with open('tfs_jaspar_ids.txt', 'w') as f:
            for tf in tfs:
                f.write(f'{tf}\\n')
    """

    stub:
    """
    touch tfs_jaspar_ids.txt
    """
}