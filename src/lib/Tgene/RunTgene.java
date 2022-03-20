package lib.Tgene;

import lib.BinaryTree.ChromosomeRegionTrees;
import lib.ExecutableStep;
import lib.GeneRegion;
import lib.Region;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class RunTgene extends ExecutableStep
{
    private Config<File> d_input;
    private final Config<File> f_output_gtf = TFPRIO.configs.tgene.fileStructure.f_transcripts_gtf;
    private final Config<File> d_output = TFPRIO.configs.tgene.fileStructure.d_output;
    private final Config<File> f_geneAnnotation = TFPRIO.configs.tepic.geneAnnotationFile;
    private final Config<File> pathToTgeneExecutable = TFPRIO.configs.tgene.pathToExecutable;

    private final Config<String> mtWriting = TFPRIO.configs.tgene.mtWriting;
    private final Config<Boolean> noClosestLocus = TFPRIO.configs.tgene.noClosestLocus;
    private final Config<Boolean> noClosestTss = TFPRIO.configs.tgene.noClosestTss;
    private final Config<Integer> maxLinkDistance = TFPRIO.configs.tgene.maxLinkDistance;
    private final Config<Double> pValue = TFPRIO.configs.tgene.pValue;

    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, f_geneAnnotation, pathToTgeneExecutable));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(f_output_gtf, d_output));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(Arrays.asList(mtWriting, noClosestLocus, noClosestTss, maxLinkDistance, pValue));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = TFPRIO.latestInputDirectory;
    }

    @Override protected void execute()
    {
        logger.info("Consensus approach with TGene is used. Preprocessing ...");

        //check_tepic_input_with_options();

        logger.info("Used data: " + d_input);

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

        ChromosomeRegionTrees trees = new ChromosomeRegionTrees();
        trees.addAllOptimized(regionsMerged);


        File tgeneExecutable = extend(pathToTgeneExecutable.get(), "bin", "tgene");

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    try
                    {
                        // Prevent crashes due to empty input files
                        if (isEmpty(f_sample))
                        {
                            continue;
                        }
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                        System.exit(1);
                    }
                    String[] name_split = f_sample.getName().split("\\.");

                    File d_output_sample = extend(d_output.get(), d_group.getName(), d_hm.getName(), name_split[0]);
                    try
                    {
                        makeSureDirectoryExists(d_output_sample);
                    } catch (IOException e)
                    {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    File gtf = f_output_gtf.get();

                    String command_execute = tgeneExecutable.getAbsolutePath();
                    command_execute += " " + f_sample.getAbsolutePath();
                    command_execute += " " + gtf.getAbsolutePath();

                    command_execute += " -oc " + d_output_sample;

                    if (noClosestLocus.get())
                    {
                        command_execute += " --no-closest-locus";
                    }
                    if (noClosestTss.get())
                    {
                        command_execute += " --no-closest-tss";
                    }

                    command_execute += " --max-link-distances " + maxLinkDistance.get();
                    command_execute += " --max-pvalue " + pValue.get();

                    //now execute TGENE:
                    String finalCommand_execute = command_execute;
                    executorService.submit(() ->
                    {
                        try
                        {
                            executeAndWait(finalCommand_execute, logger);
                        } catch (IOException | InterruptedException e)
                        {
                            e.printStackTrace();
                            System.exit(1);
                        }
                    });
                }
            }
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
