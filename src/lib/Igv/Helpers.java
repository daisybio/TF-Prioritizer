package lib.Igv;

import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;
import util.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static util.FileManagement.extend;

public class Helpers
{
    static Map<String, String> getGeneCoordinates(File f_input, Logger logger)
    {
        Map<String, String> geneCoordinates = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
        {
            String inputLine;
            reader.readLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");

                if (split[1].startsWith("CHR"))
                {
                    continue;
                }

                String gene_name;
                if (split[0].equals(""))
                {
                    gene_name = split[6];
                } else
                {
                    gene_name = split[0];
                }

                String chr = "chr" + split[1];
                int begin = Integer.parseInt(split[2]);
                int end = Integer.parseInt(split[3]);

                begin -= 50000;
                end += 50000;

                String position_string = chr + ":" + begin + "-" + end;
                geneCoordinates.put(gene_name.toUpperCase(), position_string);
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        return geneCoordinates;
    }

    static List<String> getInputFiles(Iterable<String> groups, AbstractConfig<List<String>> includePredictionData,
                                      AbstractConfig<File> d_input_tepic, AbstractConfig<File> pathToTfChipSeq,
                                      AbstractConfig<File> pathToTdf, AbstractConfig<File> d_input_peakFiles,
                                      List<File> tdfFiles)
    {
        List<String> loadFiles = new ArrayList<>();

        for (String group : groups)
        {
            //include predictive HMs
            if (includePredictionData.isSet() && includePredictionData.get().size() > 0)
            {
                for (Object includingHmObject : includePredictionData.get())
                {
                    String includingHm = (String) includingHmObject;

                    File d_input = extend(d_input_tepic.get(), group, includingHm);

                    if (d_input.exists())
                    {
                        for (File f_input : Objects.requireNonNull(d_input.listFiles(Filters.fileFilter)))
                        {
                            loadFiles.add(f_input.getAbsolutePath());
                        }
                    }
                }
            }

            //add all own TF ChIP-seq if available
            if (pathToTfChipSeq.isSet())
            {
                File d_input = extend(pathToTfChipSeq.get(), group);
                if (d_input.exists() && d_input.isDirectory())
                {
                    for (File d_input_tf : Objects.requireNonNull(d_input.listFiles(Filters.directoryFilter)))
                    {
                        for (File f_input : Objects.requireNonNull(d_input_tf.listFiles(Filters.fileFilter)))
                        {
                            loadFiles.add(f_input.getAbsolutePath());
                        }
                    }
                }
            }


            //add TDF if available
            if (pathToTdf.isSet())
            {
                File d_input = extend(pathToTdf.get(), group);
                if (d_input.exists() && d_input.isDirectory())
                {
                    for (File d_input_hm : Objects.requireNonNull(d_input.listFiles(Filters.directoryFilter)))
                    {
                        for (File f_input : Objects.requireNonNull(d_input_hm.listFiles(Filters.fileFilter)))
                        {
                            loadFiles.add(f_input.getAbsolutePath());
                            tdfFiles.add(f_input);
                        }
                    }
                }
            }
        }

        if (d_input_peakFiles.get().exists())
        {
            for (File d_input_tf : Objects.requireNonNull(d_input_peakFiles.get().listFiles(Filters.directoryFilter)))
            {
                for (File f_input : Objects.requireNonNull(d_input_tf.listFiles(Filters.fileFilter)))
                {
                    String file_name = f_input.getName();

                    if (file_name.endsWith("idx") || file_name.endsWith("bam"))
                    {
                        continue;
                    }

                    loadFiles.add(f_input.getAbsolutePath());
                }
            }
        }
        return loadFiles;
    }
}
