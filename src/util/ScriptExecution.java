package util;

import tfprio.tfprio.TFPRIO;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ScriptExecution
{
    public static void executeAndWait(File file, Logger logger)
    {
        Process process = execute(file, logger);
        waitFor(process, logger, List.of(file.getAbsolutePath()));
    }

    public static void executeAndWait(List<String> command, Logger logger)
    {
        Process process = execute(command, logger, false);
        waitFor(process, logger, command);
    }

    private static Process executeProcessBuilder(ProcessBuilder builder, Logger logger, boolean redirectOutput)
    {
        logger.debug("Executing command: " + builder.command());
        if (redirectOutput ||
                (TFPRIO.configs != null && TFPRIO.configs.general.redirectExternalScriptOutputStream.get()))
        {
            builder.redirectOutput(ProcessBuilder.Redirect.INHERIT);
        }
        if (redirectOutput ||
                (TFPRIO.configs != null && TFPRIO.configs.general.redirectExternalScriptErrorStream.get()))
        {
            builder.redirectError(ProcessBuilder.Redirect.INHERIT);
        }
        try
        {
            return builder.start();
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        // Will never be reached since logger.error exits the program.
        return null;
    }

    public static void executeAndWait(String command, Logger logger)
    {
        executeAndWait(command, logger, false);
    }

    public static void executeAndWait(String command, Logger logger, boolean redirectOutput)
    {
        Process process = execute(command, logger, redirectOutput);
        int returnCode = waitFor(process, logger, new ArrayList<>(List.of(command)));

        if (returnCode != 0)
        {
            logger.error("Received return code " + returnCode + "\n\n Command was: " + command);
        }
    }

    public static void executeAndWait(String executable, String fileExtension, Logger logger)
    {
        List<String> command = getExecutionCommand(executable, fileExtension);

        executeAndWait(command, logger);
    }

    private static int waitFor(Process process, Logger logger, List<String> command)
    {
        try
        {
            int returnCode = process.waitFor();
            if (returnCode != 0)
            {
                logger.warn("Received return code " + returnCode + " Command was: " + command);
                logger.error(new String(process.getErrorStream().readAllBytes()));
            }
            return returnCode;
        } catch (InterruptedException | IOException e)
        {
            logger.error(e.getMessage());
        }

        // Will never be reached since logger.error exits the process
        return 1;
    }

    private static List<String> getExecutionPrefix(String fileExtension, boolean fileExecution)
    {
        List<String> command = new ArrayList<>();
        switch (fileExtension)
        {
            case ".R" -> command.add("Rscript");
            case ".py" ->
            {
                command.add("python3");
                if (!fileExecution)
                {
                    command.add("-c");
                }
            }
            case ".sh" ->
            {
                command.add("sh");
            }
            default -> throw new RuntimeException("This file type is not supported.");
        }

        return command;
    }

    private static List<String> getExecutionCommand(String executable, String fileExtension)
    {
        List<String> command = getExecutionPrefix(fileExtension, false);
        command.add(executable);
        return command;
    }

    private static List<String> getExecutionCommand(File file)
    {
        String extension = file.getName().substring(file.getName().lastIndexOf("."));
        List<String> command = getExecutionPrefix(extension, true);
        command.add(file.getAbsolutePath());
        return command;
    }

    public static Process execute(String command, Logger logger)
    {
        return execute(command, logger, false);
    }

    public static Process execute(String command, Logger logger, boolean redirectOutput)
    {
        Process process = executeProcessBuilder(new ProcessBuilder(command.split(" ")), logger, redirectOutput);
        assert process != null;
        return process;
    }

    public static Process execute(List<String> command, Logger logger, boolean redirectOutput)
    {
        Process process = executeProcessBuilder(new ProcessBuilder(command), logger, redirectOutput);
        assert process != null;
        return process;
    }

    public static Process execute(List<String> command, Logger logger)
    {
        return execute(command, logger, false);
    }

    public static Process execute(File file, Logger logger)
    {
        List<String> command = getExecutionCommand(file);
        return execute(command, logger);
    }
}
