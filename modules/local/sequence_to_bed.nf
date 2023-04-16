process SEQUENCE_TO_BED {
    container 'ubuntu:22.04'

    input:
    tuple val(hm), val(group), path(binding_sequences)

    output:
    tuple val(hm), val(group), path("*.bed")

    script:
    """
    binding_sequences.sh $binding_sequences
    """
}