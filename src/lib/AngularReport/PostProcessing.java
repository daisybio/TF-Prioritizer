package lib.AngularReport;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class PostProcessing extends ExecutableStep
{
    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>();
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>();
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        boolean compress = false;

        File f_index = extend(TFPRIO.configs.angularReport.fileStructure.d_output.get(), "index.html");
        try
        {
            String content = readFile(f_index);
            writeFile(f_index, content.replace(" type=\"module\"", ""));
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        if (compress)
        {
            File d_assetsInputData =
                    extend(TFPRIO.configs.angularReport.fileStructure.d_output.get(), "assets", "input", "data");
            for (File d_tf : Objects.requireNonNull(d_assetsInputData.listFiles(
                    file -> !Arrays.asList("importantLoci", "topLog2fc").contains(file.getName()))))
            {
                File d_igv = extend(d_tf, "validation", "igv");
                compressPngsInDirectory(d_igv);
            }
        }
    }

    private void compressPngsInDirectory(File d_parent)
    {
        for (File d_child : Objects.requireNonNull(d_parent.listFiles(Filters.directoryFilter)))
        {
            compressPngsInDirectory(d_child);
        }

        for (File f_png : Objects.requireNonNull(d_parent.listFiles(Filters.getSuffixFilter(".png"))))
        {
            executorService.submit(() -> compressPng(f_png));
        }
    }

    private void compressPng(File f_png)
    {
        File f_temp = extend(f_png.getParentFile(), "temp_" + f_png.getName());

        copyFile(f_png, f_temp, logger);
        deleteFileStructure(f_png, logger);

        String command =
                "pngtopnm " + f_temp.getAbsolutePath().replace(" ", "\\ ") + " | pnmquant 16 | pnmtopng >" + " " +
                        f_png.getAbsolutePath().replace(" ", "\\ ");
        String[] cmd = {"/bin/sh", "-c", command};
        executeAndWait(List.of(cmd), logger);
        deleteFileStructure(f_temp, logger);
    }
}
