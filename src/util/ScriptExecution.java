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

        try
        {
            executeAndWait(command);
        } catch (InterruptedException | IOException e)
        {
            logger.error("Execution of command was not successful. Message: " + e.getMessage());
        }
    }

    private static void executeAndWait(List<String> command) throws IOException, InterruptedException
    {
        ProcessBuilder pb = new ProcessBuilder(command);
        if (TFPRIO.configs.general.redirectExternalScriptOutputStream.get())
        {
            pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
        }
        if (TFPRIO.configs.general.redirectExternalScriptErrorStream.get())
        {
            pb.redirectError(ProcessBuilder.Redirect.INHERIT);
        }
        Process process = pb.start();
        int returnCode = process.waitFor();

        if (returnCode != 0)
        {
            throw new ExternalScriptException(returnCode, process.getErrorStream().toString());
        }
    }

    public static void executeAndWait(String command, Logger logger) throws IOException, InterruptedException
    {
        logger.debug("Executing command: " + command);
        Process child = Runtime.getRuntime().exec(command);

        int code = child.waitFor();
        if (code != 0)
        {
            String message = child.getErrorStream().toString();
            throw new ExternalScriptException(code, message);
        }
    }

    public static void executeAndWait(String executable, String fileExtension)
    {
        List<String> command = getExecutionCommand(executable, fileExtension);
        try
        {
            executeAndWait(command);
        } catch (InterruptedException | IOException e)
        {
            System.out.println("Execution of command was not successful. Message: " + e.getMessage());
        }
    }

    private static List<String> getExecutionPrefix(String fileExtension, boolean fileExecution)
    {
        List<String> command = new ArrayList<>();
        switch (fileExtension)
        {
            case ".R" -> command.add("Rscript");
            case ".py" -> {
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
