package util.Report;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Objects;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class FileManagement
{
    static final ExecutorService executorService = Executors.newFixedThreadPool(10);

    static String loadFile(String path) throws IOException
    {
        File source = new File(path);
        return Files.readString(source.toPath());
    }

    static void copyFile(File source, File target)
    {
        copyFile(source, target, true);
    }

    static void copyFile(File source, File target, boolean compression)
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
                executorService.submit(() ->
                {
                    try
                    {
                        Process child = Runtime.getRuntime().exec(cmd);
                        child.waitFor();

                    } catch (IOException | InterruptedException e)
                    {
                        e.printStackTrace();
                    }

                });
            } else
            {
                executorService.submit(() ->
                {
                    try
                    {
                        Files.copy(source.toPath(), target.toPath(), REPLACE_EXISTING);
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                    }
                });
            }
        }
    }

    static void copyDirectory(File source, File target, boolean compression)
    {
        for (File sourceFile : Objects.requireNonNull(source.listFiles()))
        {
            File targetFile = new File(target.getAbsolutePath() + File.separator + sourceFile.getName());

            copyFile(sourceFile, targetFile, compression);
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
}
