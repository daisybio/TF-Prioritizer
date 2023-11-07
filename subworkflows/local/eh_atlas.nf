include { FETCH_LINKS } from "../../modules/local/eh_atlas/fetch_links"
include { LIFTOVER } from "../../modules/local/liftover"

def get_eh_atlas_genome(tax_id) {
    return [
        9606: "hs",
        10090: "mm",
        10116: "rn",
        7955: "dr",
        7227: "dm",
        6239: "ce",
        4932: "sc",
        9031: "gg",
        9823: "ss"
    ][tax_id]
}

def get_eh_atlas_version(species) {
    return [
        "hs": "hg19",
        "mm": "mm9",
        "dr": "danRer10",
        "dm": "dm3",
        "ce": "ce10",
        "rn": "rn5",
        "sc": "sacCer3",
        "gg": "galGal4",
        "ss": "susScr3"
    ][species]
}

workflow EH_ATLAS {
    take:
        target_genome
        tax_id

    main:
        eh_atlas_genome = get_eh_atlas_genome(tax_id)
        eh_atlas_version = get_eh_atlas_version(eh_atlas_genome)

        FETCH_LINKS(eh_atlas_genome)

        // LIFTOVER(ch_eh_atlas, eh_atlas_version, target_genome)

    //emit:
        //bed = LIFTOVER.out
}
