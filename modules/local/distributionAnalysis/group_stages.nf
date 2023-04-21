process GROUP_STAGES {
    conda 'pandas'
    container 'tfprio-python'

    input:
        tuple val(hm), val(same), path(files)

    output:
        tuple val(hm), val(same), path("${hm}_${same}.tsv")

    script:
        """
        join.py --input ${files} --output ${hm}_${same}.tsv
        """
}