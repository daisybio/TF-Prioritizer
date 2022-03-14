package util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ScriptExecution
{
    public static int executeAndWait(File file, Logger logger) throws InterruptedException, IOException
    {
        List<String> command = getExecutionCommand(file);
        logger.info("Executing command: " + command);

        return executeAndWait(command);
    }

    private static int executeAndWait(List<String> command) throws IOException, InterruptedException
    {
        ProcessBuilder pb = new ProcessBuilder(command);
        pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
        pb.redirectError(ProcessBuilder.Redirect.INHERIT);
        Process process = pb.start();
        return process.waitFor();
    }

    public static int executeAndWait(String executable, String fileExtension)
    {
        List<String> command = getExecutionCommand(executable, fileExtension);
        int ret;
        try
        {
            ret = executeAndWait(command);
        } catch (IOException | InterruptedException e)
        {
            ret = 1;
        }
        return ret;
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
