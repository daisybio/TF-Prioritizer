include { DOWNLOAD } from "../../modules/local/chip_atlas/download"

workflow CHIP_ATLAS {
    main:
        DOWNLOAD()

}