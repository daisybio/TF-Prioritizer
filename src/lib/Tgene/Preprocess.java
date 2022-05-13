package lib.Tgene;

import lib.ExecutableStep;
import util.Regions.GeneRegion;
import util.Regions.Region;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.*;
import java.util.*;

import static util.FileManagement.makeSureFileExists;

public class Preprocess extends ExecutableStep
{
    private final AbstractConfig<File> f_geneAnnotation = TFPRIO.configs.tepic.geneAnnotationFile;

    private final GeneratedFileStructure f_output_gtf = TFPRIO.configs.tgene.fileStructure.f_transcripts_gtf;
    private final GeneratedFileStructure f_output_regions = TFPRIO.configs.tgene.fileStructure.f_preprocessing_regions;

    private final AbstractConfig<String> mtWriting = TFPRIO.configs.tgene.mtWriting;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(List.of(f_geneAnnotation));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_output_gtf, f_output_regions));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(mtWriting));
    }

    @Override protected void execute()
    {
        File f_output = f_output_gtf.get();
        try
        {
            makeSureFileExists(f_output);
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        List<GeneRegion> regions = new ArrayList<>();

        try (BufferedReader transcriptsReader = new BufferedReader(new FileReader(f_geneAnnotation.get()));
             BufferedWriter transcriptsWriter = new BufferedWriter(new FileWriter(f_output)))
        {
            String inputLine;

            while ((inputLine = transcriptsReader.readLine()) != null)
            {
                if (inputLine.startsWith("#"))
                {
                    continue;
                }

                String[] split = inputLine.split("\t");

                String type = split[2];
                String chromosome = split[0];

                if (type.equals("transcript"))
                {
                    if (chromosome.matches(".*M.*"))
                    {
                        chromosome = mtWriting.get();
                    }
                    if (chromosome.startsWith("chr"))
                    {
                        chromosome = chromosome.substring("chr".length());
                    }

                    String[] contentSplitted = split[8].split(";");
                    String geneID = "";
                    for (String entry : contentSplitted)
                    {
                        if (entry.startsWith("gene_id"))
                        {
                            geneID = entry.substring(entry.indexOf("\"") + 1, entry.lastIndexOf("\""));
                            geneID = geneID.substring(0, geneID.indexOf("."));
                            break;
                        }
                    }

                    regions.add(
                            new GeneRegion(chromosome, Integer.parseInt(split[3]), Integer.parseInt(split[4]), geneID));

                    StringBuilder sb_line = new StringBuilder(chromosome);
                    for (int i = 1; i < split.length; i++)
                    {
                        sb_line.append("\t").append(split[i]);
                    }
                    transcriptsWriter.write(sb_line.toString());
                    transcriptsWriter.newLine();
                }
            }
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }
        List<Region> regionsMerged = merge(regions);

        makeSureFileExists(f_output_regions.get(), logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_regions.get())))
        {
            for (Region region : regionsMerged)
            {
                writer.write(region.toString());
                writer.newLine();
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    private List<Region> merge(List<GeneRegion> inputRegions)
    {
        List<Region> outputRegions = new ArrayList<>();

        for (int i = 0; i < inputRegions.size(); i++)
        {
            GeneRegion region = inputRegions.get(i);
            for (int j = i + 1; j < inputRegions.size(); j++)
            {
                GeneRegion other = inputRegions.get(j);

                if (region.getId().equals(other.getId()))
                {
                    region.merge(other);
                } else
                {
                    i = j - 1;
                    break;
                }
            }
            outputRegions.add(region);
        }
        return outputRegions;
    }
}
