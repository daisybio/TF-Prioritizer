package util;

import tfprio.TFPRIO;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Logger
{
    private final SimpleDateFormat formatter = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");

    boolean fileLoggingEnabled;
    File logFile;
    String module;

    private enum LogLevel
    {
        DEBUG, INFO, WARN, ERROR, PROG
    }

    public Logger(String module)
    {
        this(module, TFPRIO.configs.general.fileLogging.get(), TFPRIO.configs.general.logFile.get());
    }

    public Logger(String module, boolean fileLoggingEnabled, File logFile)
    {
        this.fileLoggingEnabled = fileLoggingEnabled;
        this.module = module;

        if (fileLoggingEnabled)
        {
            assert logFile != null;
            this.logFile = logFile;
        }

        debug("Initialised logger " +
                (fileLoggingEnabled ? "with file logging to " + logFile.getAbsolutePath() : "without file logging."));
    }

    public void debug(String message)
    {
        logLine(message, LogLevel.DEBUG);
    }

    public void info(String message)
    {
        logLine(message, LogLevel.INFO);
    }

    public void warn(String message)
    {
        logLine(message, LogLevel.WARN);
    }

    public void error(String message)
    {
        logLine(message, LogLevel.ERROR);
        System.exit(1);
    }

    private String getLine(String message, LogLevel level)
    {
        return "[" + formatter.format(new Date()) + "]\t" + level + "\t[" + module + "]\t" + message;
    }

    private void logLine(String message, LogLevel level)
    {
        String line = getLine(message, level);

        if (fileLoggingEnabled)
        {
            try
            {
                FileManagement.appendToFile(logFile, line + "\n");
            } catch (IOException e)
            {
                fileLoggingEnabled = false;
                logLine("File logging not possible: " + e.getMessage(), LogLevel.ERROR);
            }
        }

        System.out.println(line);
    }

    public void progress(double percentage)
    {
        String message = String.format("Progress: %.2f %%", percentage);
        System.out.print(getLine(message, LogLevel.PROG) + "\r");
    }

    @Deprecated public void logLine(String message)
    {
        logLine(message, LogLevel.INFO);
    }
}
