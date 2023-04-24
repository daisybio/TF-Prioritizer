process REPORT {
    input:
        path(data)
    
    output:
        path("report")

    script:
        """
        cp $projectDir/assets/report report
        """
}