process REPORT {
    container 'tfprio-angular'

    input:
        path(data)
    
    output:
        path("tfprio_report")

    script:
        """
        cp -r $projectDir/assets/report report
        cd report

        npm install
        ng build --output-path ../tfprio_report

        cp $projectDir/bin/serve_report.sh ../tfprio_report
        """
}