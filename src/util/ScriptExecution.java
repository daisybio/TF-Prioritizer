package util;

import java.io.File;
import java.io.IOException;

public class ScriptExecution
{
    public static int executeAndWait(File file, Logger logger) throws InterruptedException, IOException
    {
        String command = getExecutionCommand(file);
        logger.info("Executing command: " + command);

        Process process = Runtime.getRuntime().exec(command);
        return process.waitFor();
    }

    private static String getExecutionCommand(File file)
    {
        String extension = file.getName().substring(file.getName().lastIndexOf("."));
        String prefix = switch (extension)
                {
                    case ".R" -> "Rscript ";
                    case ".py" -> "python3 ";
                    default -> throw new RuntimeException("This file type is not supported.");
                };

        return prefix + file.getAbsolutePath();
    }
}
