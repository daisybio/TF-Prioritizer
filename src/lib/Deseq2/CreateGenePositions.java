package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreateGenePositions extends ExecutableStep
{
    private final AbstractConfig<File> d_input_enhancerDbs = TFPRIO.configs.deSeq2.d_enhancerDB;
    private final GeneratedFileStructure f_script =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_script;
    private final GeneratedFileStructure f_data_prev =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_genePositionsPrev;
    private final GeneratedFileStructure f_data_version =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_version;
    private final GeneratedFileStructure d_output_enhancerDbs =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_genePositions_enhancerDBs;
    private final GeneratedFileStructure f_output_mergedEnhancerDbs =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_mergedEnhancerDbs;
    private final AbstractConfig<File> f_mapping = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final GeneratedFileStructure f_data =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_data;
    private final AbstractConfig<File> f_upliftScriptTemplate =
            TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingUplift;
    private final GeneratedFileStructure f_upliftScript =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_uplift;

    private final AbstractConfig<File> f_scriptTemplate =
            TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingGetGenePositions;
    private final AbstractConfig<File> f_dbScriptTemplate =
            TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingDbUplift;

    private final AbstractConfig<Boolean> calculateGenePositionsEnabled =
            TFPRIO.configs.general.calculateGenePositionsEnabled;
    private final AbstractConfig<String> speciesReferenceGenome = TFPRIO.configs.igv.speciesReferenceGenome;
    private final AbstractConfig<String> speciesBiomart = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final AbstractConfig<Map<String, String>> synonymDict = TFPRIO.configs.igv.grcSynonymDict;
    private final AbstractConfig<String> bedFormat =
            TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedFormat;
    private final AbstractConfig<String> bedSuffix =
            TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix;

    private final AbstractConfig<List<String>> enhancerDBs = TFPRIO.configs.igv.enhancerDatabases;


    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_mapping, f_upliftScriptTemplate))
        {{
            if (enhancerDBs.isSet())
            {
                add(d_input_enhancerDbs);
                add(f_dbScriptTemplate);
            }
        }};
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_script, f_data_prev, f_data_version, f_data, f_upliftScript))
        {{
            for (GeneratedFileStructure notGeneratedFileStructure : Arrays.asList(d_output_enhancerDbs,
                    f_output_mergedEnhancerDbs))
            {
                if (enhancerDBs.isSet())
                {
                    add(notGeneratedFileStructure);
                } else
                {
                    notGeneratedFileStructure.deleteAndSetNoGenerationReason(enhancerDBs.getName() + " not set");
                }
            }
        }};
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(calculateGenePositionsEnabled, speciesBiomart, speciesReferenceGenome, synonymDict))
        {{
            if (enhancerDBs.isSet())
            {
                add(bedFormat);
                add(bedSuffix);
            }
        }};
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(List.of(enhancerDBs));
    }

    @Override protected void execute()
    {
        //get biomart gene positions
        if (calculateGenePositionsEnabled.get())
        {
            try
            {
                String script = readFile(f_scriptTemplate.get());
                script = script.replace("{INPUTFILE}", f_mapping.get().getAbsolutePath());
                script = script.replace("{SPECIES}", speciesBiomart.get());
                script = script.replace("{DATA_PREV_FILE}", f_data_prev.get().getAbsolutePath());
                script = script.replace("{VERSIONFILE}", f_data_version.get().getAbsolutePath());
                writeFile(f_script.get(), script);
            } catch (IOException e)
            {
                e.printStackTrace();
            }
            executeAndWait(f_script.get(), logger);
        } else
        {
            logger.info("Gene positions were not fetched since command line parameter -b was set.");
        }

        try
        {
            makeSureFileExists(f_data.get());

            String script = readFile(f_upliftScriptTemplate.get());
            script = script.replace("{INPUTFILE}", f_data_prev.get().getAbsolutePath());
            script = script.replace("{VERSIONFILE}", f_data_version.get().getAbsolutePath());
            script = script.replace("{OUTPUTFILE}", f_data.get().getAbsolutePath());
            script = script.replace("{REFERENCE_GENOME}", speciesReferenceGenome.get());

            StringBuilder sb_grc = new StringBuilder();

            for (String key : synonymDict.get().keySet())
            {
                sb_grc.append("'").append(key).append("':");
                sb_grc.append("'").append(synonymDict.get().get(key)).append("',");
            }
            sb_grc.deleteCharAt(sb_grc.lastIndexOf(","));

            script = script.replace("{REFERENCE_GENOME_DICT}", sb_grc.toString());
            writeFile(f_upliftScript.get(), script);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        String line_biomart_version = null;
        try (BufferedReader versionReader = new BufferedReader(new FileReader(f_data_version.get())))
        {
            line_biomart_version = versionReader.readLine();
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }


        String biomartVersion = line_biomart_version.split("\\.")[0];

        logger.info("Uplift positions to correct genome version");
        logger.info("Our genome: " + speciesBiomart.get());

        if (synonymDict.get().containsKey(biomartVersion))
        {
            logger.info("Genome got from biomaRt: " + biomartVersion + " = " + synonymDict.get().get(biomartVersion));

            if (speciesBiomart.get().equals(synonymDict.get().get(biomartVersion)))
            {
                logger.info("equal genome versions. No uplifted needed. Skipping this step.");
            } else
            {
                try
                {
                    Thread.sleep(1000);
                } catch (InterruptedException ignore)
                {
                }
                executeAndWait(f_upliftScript.get(), logger);
            }
        } else
        {
            logger.error("Do not have " + biomartVersion + " in memory. Please add at igv_GRC_synonym_dict");
        }

        if (enhancerDBs.isSet())
        {
            String scriptTemplate = null;
            try
            {
                scriptTemplate = readFile(f_dbScriptTemplate.get());
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
            assert scriptTemplate != null;

            String speciesSymbol = speciesReferenceGenome.get().replaceAll("[^A-Za-z]", "");

            StringBuilder sb_merged = new StringBuilder();
            sb_merged.append(bedFormat.get()).append("\n");

            for (String enhancerDB : enhancerDBs.get())
            {
                File f_input = extend(d_input_enhancerDbs.get(), enhancerDB + bedSuffix.get());
                File f_output = extend(d_output_enhancerDbs.get(), enhancerDB + "_prev" + bedSuffix.get());

                StringBuilder sb_specific = new StringBuilder();
                sb_specific.append(bedFormat.get()).append("\n");

                boolean needUplift = true;
                String lastFoundReferenceGenome = "";

                try (BufferedReader reader = new BufferedReader(new FileReader(f_input)))
                {
                    String inputLine;
                    reader.readLine();

                    while ((inputLine = reader.readLine()) != null)
                    {
                        String[] split = inputLine.split("\t");
                        String foundReferenceGenome = split[4];
                        String foundSpeciesSymbol = foundReferenceGenome.replaceAll("[^A-Za-z]", "");

                        if (foundSpeciesSymbol.equals(speciesSymbol))
                        {
                            if (foundReferenceGenome.equals(speciesReferenceGenome.get()))
                            {
                                sb_merged.append(inputLine);
                                sb_merged.append("\n");
                                needUplift = false;
                                continue;
                            }

                            if (lastFoundReferenceGenome.equals(""))
                            {
                                lastFoundReferenceGenome = foundReferenceGenome;
                            }

                            if (!foundReferenceGenome.equals(lastFoundReferenceGenome))
                            {
                                continue;
                            }

                            sb_specific.append(inputLine);
                            sb_specific.append("\n");
                        }
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }

                if (needUplift)
                {
                    try
                    {
                        writeFile(f_output, sb_specific.toString());
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    File f_upliftedDB = extend(d_output_enhancerDbs.get(), enhancerDB + bedSuffix.get());

                    File f_dbScript = extend(d_output_enhancerDbs.get(), enhancerDB + "_" +
                            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_uplift.get().getName());

                    String script = scriptTemplate;

                    script = script.replace("{INPUTFILE}", f_input.getAbsolutePath());
                    script = script.replace("{OUTPUTFILE}", f_upliftedDB.getAbsolutePath());
                    script = script.replace("{REFERENCE_GENOME}", speciesReferenceGenome.get());
                    script = script.replace("{FOUND_GENOME}", lastFoundReferenceGenome);

                    try
                    {
                        writeFile(f_dbScript, script);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    logger.info("Uplift db " + enhancerDB + " positions to correct genome version");

                    executeAndWait(f_dbScript, logger);

                    try (BufferedReader reader = new BufferedReader(new FileReader(f_upliftedDB.getAbsolutePath())))
                    {
                        String inputLine;
                        reader.readLine();
                        while ((inputLine = reader.readLine()) != null)
                        {
                            sb_merged.append(inputLine);
                            sb_merged.append("\n");
                        }
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }
                }
            }

            try
            {
                writeFile(f_output_mergedEnhancerDbs.get(), sb_merged.toString());
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
        }
    }
}
