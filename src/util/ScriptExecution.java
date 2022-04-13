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
        logger.debug("Executing command: " + command);

        executeAndWait(command, logger);
    }

    public static void executeAndWait(List<String> command, Logger logger)
    {
        logger.debug("Executing command: " + command);
        ProcessBuilder pb = new ProcessBuilder(command);
        if (TFPRIO.configs.general.redirectExternalScriptOutputStream.get())
        {
            pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
        }
        if (TFPRIO.configs.general.redirectExternalScriptErrorStream.get())
        {
            pb.redirectError(ProcessBuilder.Redirect.INHERIT);
        }
        try
        {
            Process process = pb.start();
            int returnCode = process.waitFor();

            if (returnCode != 0)
            {
                String message = new String(process.getErrorStream().readAllBytes());
                logger.error("Received return code " + returnCode + ":\n\n" + message);
            }
        } catch (IOException | InterruptedException e)
        {
            e.printStackTrace();
        }


    }

    public static void executeAndWait(String command, Logger logger)
    {
        logger.debug("Executing command: " + command);
        try
        {
            Process child = Runtime.getRuntime().exec(command);

            int code = child.waitFor();
            if (code != 0)
            {
                String message = new String(child.getErrorStream().readAllBytes());
                logger.error("Received return code " + code + ":\n\n" + message);
            }
        } catch (IOException | InterruptedException e)
        {
            logger.error(e.getMessage());
        }
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
