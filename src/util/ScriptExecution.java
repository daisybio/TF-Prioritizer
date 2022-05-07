package util;

import tfprio.TFPRIO;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ScriptExecution
{
    public static void executeAndWait(File file, Logger logger)
    {
        List<String> command = getExecutionCommand(file);

        executeAndWait(command, logger);
    }

    public static void executeAndWait(List<String> command, Logger logger)
    {
        ProcessBuilder pb = new ProcessBuilder(command);
        executeProcessBuilder(pb, logger, false);
    }

    private static void executeProcessBuilder(ProcessBuilder builder, Logger logger, boolean redirectOutput)
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
            Process process = builder.start();
            int returnCode = process.waitFor();

            if (returnCode != 0)
            {
                String message = new String(process.getErrorStream().readAllBytes());
                logger.error("Received return code " + returnCode + "\n\n Command was: " + builder.command() + "\n\n" +
                        message);
            }
        } catch (IOException | InterruptedException e)
        {
            logger.error(e.getMessage());
        }
    }

    public static void executeAndWait(String command, Logger logger)
    {
        executeAndWait(command, logger, false);
    }

    public static void executeAndWait(String command, Logger logger, boolean redirectOutput)
    {
        executeProcessBuilder(new ProcessBuilder(command.split(" ")), logger, redirectOutput);
    }

    public static void executeAndWait(String executable, String fileExtension, Logger logger)
    {
        List<String> command = getExecutionCommand(executable, fileExtension);

        executeAndWait(command, logger);
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
}
