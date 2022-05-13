package lib.Plots;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class OpenRegionsViolinPlots extends ExecutableStep
{
    private final AbstractConfig<File> d_input = TFPRIO.configs.tepic.fileStructure.d_outputRaw;
    private final AbstractConfig<File> f_scriptTemplate =
            TFPRIO.configs.scriptTemplates.f_plots_openChromatinViolinPlots;

    private final GeneratedFileStructure f_output_data =
            TFPRIO.configs.tepic.fileStructure.f_postprocessing_openChromatinViolins_data_csv;
    private final GeneratedFileStructure f_output_script =
            TFPRIO.configs.tepic.fileStructure.f_postprocessing_openChromatinViolins_script_R;
    private final GeneratedFileStructure f_output_plot =
            TFPRIO.configs.tepic.fileStructure.f_postprocessing_openChromatinViolins_plots_image;

    private final AbstractConfig<String> s_regionsToTargetGenes =
            TFPRIO.configs.tepic.fileStructure.s_outputRaw_regionsToTargetGenes;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_input, f_scriptTemplate));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_output_data, f_output_script, f_output_plot));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(List.of(s_regionsToTargetGenes));
    }

    @Override protected void execute()
    {
        HashMap<String, ArrayList<Integer>> hm_to_openChrLengths = new HashMap<>();

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                String hm = d_hm.getName();
                ArrayList<Integer> hm_lengths;

                if (hm_to_openChrLengths.containsKey(hm))
                {
                    hm_lengths = hm_to_openChrLengths.get(hm);
                } else
                {
                    hm_lengths = new ArrayList<>();
                    hm_to_openChrLengths.put(hm, hm_lengths);
                }

                for (File d_sample : Objects.requireNonNull(d_hm.listFiles(Filters.directoryFilter)))
                {
                    for (File f_regionsToTargets : Objects.requireNonNull(d_sample.listFiles(Filters.fileFilter)))
                    {
                        if (f_regionsToTargets.getName().equals(s_regionsToTargetGenes.get()))
                        {
                            try (BufferedReader reader = new BufferedReader(new FileReader(f_regionsToTargets)))
                            {
                                String line_regions_to_targets;
                                reader.readLine();

                                while ((line_regions_to_targets = reader.readLine()) != null)
                                {
                                    String[] split = line_regions_to_targets.split("\t");
                                    String[] regionSplit = split[0].split(":");
                                    String[] coordinatesSplit = regionSplit[1].split("-");

                                    int start = Integer.parseInt(coordinatesSplit[0]);
                                    int end = Integer.parseInt(coordinatesSplit[1]);

                                    int length = end - start;

                                    hm_lengths.add(length);
                                }
                            } catch (IOException e)
                            {
                                logger.error(e.getMessage());
                            }
                        }
                    }
                }
            }
        }


        makeSureFileExists(f_output_data.get(), logger);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(f_output_data.get())))
        {
            writer.write("HM\tLENGTH\n");
            for (String hm : hm_to_openChrLengths.keySet())
            {
                for (Integer length : hm_to_openChrLengths.get(hm))
                {
                    writer.write(hm + "\t" + length + "\n");
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        makeSureFileExists(f_output_plot.get(), logger);

        //create Rscript and execute Rscript
        try
        {
            String script = readFile(f_scriptTemplate.get());
            script = script.replace("{INPUTFILE}", f_output_data.get().getAbsolutePath());
            script = script.replace("{OUTPUTFILE}", f_output_plot.get().getAbsolutePath());
            writeFile(f_output_script.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        executorService.submit(() -> executeAndWait(f_output_script.get(), logger));
    }
}
