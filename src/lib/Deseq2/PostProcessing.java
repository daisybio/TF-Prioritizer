package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;

public class PostProcessing extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.deSeq2.fileStructure.d_outputRaw;
    private final GeneratedFileStructure d_output = TFPRIO.configs.deSeq2.fileStructure.d_output;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(d_input));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>();
    }

    @Override protected void execute()
    {
        for (File f_input : Objects.requireNonNull(d_input.get().listFiles(Filters.fileFilter)))
        {
            File f_output = extend(d_output.get(), f_input.getName());
            try
            {
                makeSureFileExists(f_output);
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }

            try (BufferedReader reader = new BufferedReader(new FileReader(f_input));
                 BufferedWriter writer = new BufferedWriter(new FileWriter(f_output)))
            {
                writer.write("geneID\tlog2fc");
                writer.newLine();

                String inputLine;
                reader.readLine();

                while ((inputLine = reader.readLine()) != null)
                {
                    String[] split = inputLine.split("\t");

                    if (!split[2].equals("NA"))
                    {
                        String line = split[0].substring(1, split[0].length() - 1) + "\t" + split[2];

                        writer.write(line);
                        writer.newLine();
                    }
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        }
    }
}
