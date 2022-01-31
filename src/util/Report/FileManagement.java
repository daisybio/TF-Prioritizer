package util.Report;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class FileManagement
{
    static String loadFile(String path) throws IOException
    {
        File source = new File(path);
        return Files.readString(source.toPath());
    }

    static void copyFile(File source, File target) throws IOException
    {
        copyFile(source, target, true);
    }

    static void copyFile(File source, File target, boolean compression) throws IOException
    {
        if (target.exists())
        {
            return;
        }

        if (source.exists())
        {
            if (!target.getParentFile().exists())
            {
                target.getParentFile().mkdirs();
            }

            if (target.getName().endsWith(".png") && compression)
            {
                String command = "pngtopnm " + source.getAbsolutePath() + " | pnmquant 16 | pnmtopng > " +
                        target.getAbsolutePath();
                String[] cmd = {"/bin/sh", "-c", command};
                Process child = Runtime.getRuntime().exec(cmd);
                try
                {
                    child.waitFor(5, TimeUnit.SECONDS);
                } catch (InterruptedException ignored)
                {
                }
            } else
            {
                Files.copy(source.toPath(), target.toPath(), REPLACE_EXISTING);
            }
        }
    }

    static void copyDirectory(File source, File target, boolean compression) throws IOException
    {
        copyDirectory(source, target, compression, ".*", new ArrayList<>());
    }

    static void copyDirectory(File source, File target, boolean compression, String fileNameRegex,
                              List<String> removables) throws IOException
    {
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        for (File sourceFile : Objects.requireNonNull(source.listFiles()))
        {
            if (sourceFile.isFile())
            {
                if (sourceFile.getName().matches(fileNameRegex))
                {
                    String cleanName = sourceFile.getName();
                    String fileExtension = cleanName.substring(cleanName.lastIndexOf("."));

                    cleanName = cleanName.replace(fileExtension, "");

                    for (String removable : removables)
                    {
                        cleanName = cleanName.replace(removable, "");
                    }

                    while (cleanName.contains("__"))
                    {
                        cleanName = cleanName.replace("__", "_");
                    }
                    while (cleanName.startsWith("_"))
                    {
                        cleanName = cleanName.substring(1);
                    }
                    while (cleanName.endsWith("_"))
                    {
                        cleanName = cleanName.substring(0, cleanName.length() - 1);
                    }
                    cleanName = cleanName + fileExtension;
                    File targetFile = new File(target.getAbsolutePath() + File.separator + cleanName);
                    executorService.execute(() ->
                    {
                        try
                        {
                            copyFile(sourceFile, targetFile, compression);
                        } catch (IOException ignored)
                        {
                        }
                    });
                }
            } else
            {
                File subDir = new File(target.getAbsolutePath() + File.separator + sourceFile.getName());
                copyDirectory(sourceFile, subDir, compression, fileNameRegex, removables);
            }
        }
        try
        {
            executorService.shutdown();
            executorService.awaitTermination(30, TimeUnit.SECONDS);
        } catch (InterruptedException e)
        {
            System.out.println("[REPORT] Error during copy process: " + source.getAbsolutePath() + " to " +
                    target.getAbsolutePath());
        }
    }

    static void writeHTML(String path, String content, int relativationDepth) throws IOException
    {
        content = content.replace("{RELATIVATION}", (".." + File.separator).repeat(relativationDepth));

        writeFile(path, content);
    }

    static void writeFile(String path, String content) throws IOException
    {
        File target = new File(path);
        if (!target.exists())
        {
            target.getParentFile().mkdirs();
            target.createNewFile();
        }
        Files.writeString(target.toPath(), content);
    }

    static String findValueInTable(String term, int searchIndex, int resultIndex, File file, String sep,
                                   boolean ignoreCase) throws FileNotFoundException, NoSuchFieldException
    {
        if (file.isFile())
        {
            try (Scanner scanner = new Scanner(file))
            {
                while (scanner.hasNextLine())
                {
                    String line = scanner.nextLine();
                    if (line.split(sep).length > searchIndex && line.split(sep).length > resultIndex)
                    {

                        if (term.equals(line.split(sep)[searchIndex]) ||
                                ignoreCase && term.equalsIgnoreCase(line.split(sep)[searchIndex]))
                        {
                            return line.split(sep)[resultIndex];
                        }
                    }
                }
            }
        }
        throw new NoSuchFieldException();
    }

    static File getFileIfInDirectory(File directory, String fileNameRegex, boolean lookingForFiles)
    {
        if (directory == null)
        {
            return null;
        }

        for (File entry : Objects.requireNonNull(directory.listFiles()))
        {
            if (lookingForFiles != entry.isFile())
            {
                continue;
            }
            if (entry.getName().matches(fileNameRegex))
            {
                return entry;
            }
        }
        return null;
    }
}
