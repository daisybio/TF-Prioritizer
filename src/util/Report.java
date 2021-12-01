package util;

public class Report
{
    private final Logger logger;
    private final Options_intern options_intern;

    public Report(Options_intern options_intern) throws Exception
    {
        this.options_intern = options_intern;
        logger = new Logger(true, options_intern.com2pose_working_directory);
    }

    public void generate() throws Exception
    {
        logger.logLine("[REPORT] Start generating report");
        generateHome();

        logger.logLine("[REPORT] Finished generating report");
    }

    private void generateHome() throws Exception
    {
        logger.logLine("[REPORT] Start generating report home page");
        logger.logLine(options_intern.com2pose_working_directory);

        logger.logLine("[REPORT] Finished generating report home page");
    }
}
