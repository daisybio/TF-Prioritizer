package util;

import tfprio.tfprio.TFPRIO;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

import static util.FileManagement.appendToFile;

public class Logger
{
    private final SimpleDateFormat formatter = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");

    boolean fileLoggingEnabled;
    File logFile;
    String module;

    private enum LogLevel
    {
        DEBUG(TerminalColors.WHITE), INFO(TerminalColors.DEFAULT), WARN(TerminalColors.YELLOW),
        ERROR(TerminalColors.RED), PROG(TerminalColors.BLUE);

        public final String color;

        LogLevel(TerminalColors color)
        {
            this.color = color.code;
        }
    }

    private enum TerminalColors
    {
        RED("\u001B[31m"), GREEN("\u001B[32m"), YELLOW("\u001B[33m"), BLUE("\u001B[34m"), PURPLE("\u001B[35m"),
        CYAN("\u001B[36m"), WHITE("\u001B[37m"), DEFAULT("\u001B[0m");

        public final String code;

        TerminalColors(String code)
        {
            this.code = code;
        }
    }

    /**
     * Create a new Logger instance
     *
     * @param module the logger module name, is logged with each message
     */
    public Logger(String module)
    {
        this(module, TFPRIO.configs.general.fileLogging.get(), TFPRIO.configs.getLogFile());
    }

    /**
     * Alternative constructor which has to be used if the {@link TFPRIO} configs are not yet ready to use
     *
     * @param module             the logger module name, is logged with each message
     * @param fileLoggingEnabled defined if file logging should take place
     * @param logFile            the log file
     */
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

        if (fileLoggingEnabled && this.logFile.exists() && this.module.equals("TFPRIO"))
        {
            try
            {
                appendToFile(this.logFile, "###########################################################\n");
                appendToFile(this.logFile, "######################### NEW RUN #########################\n");
                appendToFile(this.logFile, "###########################################################\n");
            } catch (IOException e)
            {
                this.fileLoggingEnabled = false;
                warn("File logging not possible: " + e.getMessage());
            }
        }
    }

    /**
     * Log a debug message
     * <p>
     * Debug messages should only contain information that offers no additional value to the user.
     * Debug messages should help developers with locating bugs.
     *
     * @param message the log message
     */
    public void debug(String message)
    {
        logLine(message, LogLevel.DEBUG);
    }

    /**
     * Log an info message
     * <p>
     * Info messages should inform the user about the pipeline progress and the running processes during normal
     * execution.
     *
     * @param message the log message
     */
    public void info(String message)
    {
        logLine(message, LogLevel.INFO);
    }

    /**
     * Log a warning message
     * <p>
     * Warn messages should inform the user about irregularities during the pipeline execution.
     *
     * @param message the log message
     */
    public void warn(String message)
    {
        logLine(message, LogLevel.WARN);
    }

    /**
     * Log an error message
     * <p>
     * Error messages should inform the user about fatal errors during pipeline execution.
     * Makes the pipeline process terminate.
     *
     * @param message the log message
     */
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
                appendToFile(logFile, line + "\n");
            } catch (IOException e)
            {
                fileLoggingEnabled = false;
                warn("File logging not possible: " + e.getMessage());
            }
        }

        System.out.println(level.color + line + TerminalColors.DEFAULT.code);
    }

    /**
     * Log a progress message.
     * <p>
     * Progress messages should inform the user about the progress within a certain {@link lib.ExecutableStep} and
     * give an estimation of the remaining time.
     * Progress messages are not written to the log file.
     *
     * @param message the progress message
     */
    public void progress(String message)
    {
        //System.out.print(LogLevel.PROG.color + getLine(message, LogLevel.PROG) + TerminalColors.DEFAULT.code + "\r");
    }
}
