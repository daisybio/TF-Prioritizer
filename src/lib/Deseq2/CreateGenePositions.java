package lib.Deseq2;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import util.Configs.Config;
import util.ExternalScriptException;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class CreateGenePositions extends ExecutableStep
{
    private final Config<File> f_script = TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_script;
    private final Config<File> f_data_prev =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_genePositionsPrev;
    private final Config<File> f_data_version =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_version;
    private final Config<File> f_mapping = TFPRIO.configs.deSeq2.fileStructure.f_mapping;
    private final Config<File> f_data = TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_data;
    private final Config<File> f_upliftScriptTemplate = TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingUplift;
    private final Config<File> f_upliftScript =
            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_uplift;

    private final Config<File> f_scriptTemplate = TFPRIO.configs.scriptTemplates.f_deseq2PreprocessingGetGenePositions;

    private final Config<Boolean> calculateGenePositionsEnabled = TFPRIO.configs.general.calculateGenePositionsEnabled;
    private final Config<String> speciesReferenceGenome = TFPRIO.configs.igv.speciesReferenceGenome;
    private final Config<String> speciesBiomart = TFPRIO.configs.deSeq2.biomartDatasetSpecies;
    private final Config<Map> synonymDict = TFPRIO.configs.igv.grcSynonymDict;


    @Override protected Set<Config<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_scriptTemplate, f_mapping, f_upliftScriptTemplate));
    }

    @Override protected Set<Config<File>> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_script, f_data_prev, f_data_version, f_data, f_upliftScript));
    }

    @Override protected Set<Config<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(calculateGenePositionsEnabled, speciesBiomart, speciesReferenceGenome, synonymDict));
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
            try
            {
                executeAndWait(f_script.get(), logger);
            } catch (ExternalScriptException e)
            {
                e.printStackTrace();
            }
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

            for (Object keyObject : synonymDict.get().keySet())
            {
                String key = (String) keyObject;
                sb_grc.append("'").append(key).append("':");
                sb_grc.append("'").append(synonymDict.get().get(key)).append("',");
            }
            sb_grc.deleteCharAt(sb_grc.lastIndexOf(","));

            script = script.replace("{REFERENCE_GENOME_DICT}", sb_grc.toString());
            writeFile(f_upliftScript.get(), script);
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        String line_biomart_version = null;
        try (BufferedReader versionReader = new BufferedReader(new FileReader(f_data_version.get())))
        {
            line_biomart_version = versionReader.readLine();
        } catch (IOException e)
        {
            e.printStackTrace();
            System.exit(1);
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
                    executeAndWait(f_upliftScript.get(), logger);
                } catch (ExternalScriptException e)
                {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        } else
        {
            logger.info("Do not have " + biomartVersion + " in memory. Please add at igv_GRC_synonym_dict");
            System.exit(1);
        }

        // TODO: There is no data available for testing the following section. This is why i left it in the old style
        //  for now.
        /*
        if (TFPRIO.configs.igv.enhancerDatabases.get().size() > 0)
        {
            logger.info("Check reference genomes of Enhances DBs to currently used genomes and uplift if necessary.");

            String igv_species_symbol = TFPRIO.configs.igv.speciesReferenceGenome.get().replaceAll("[^A-Za-z]", "");


            File f_output_enhancer_dbs =
                    TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_genePositions_enhancerDBs.get();

            File f_input_enhancerdbs_root = TFPRIO.configs.deSeq2.d_enhancerDB.get();

            StringBuilder sb_merged_dbs = new StringBuilder();
            sb_merged_dbs.append(
                            TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedFormat.get())
                    .append("\n");

            for (Object enhancerDbObject : TFPRIO.configs.igv.enhancerDatabases.get())
            {
                String enhancer_db = (String) enhancerDbObject;

                File f_input_enhancerdbs_DB = extend(f_input_enhancerdbs_root, enhancer_db +
                        TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());

                File f_output_filtered_db = extend(f_output_enhancer_dbs, enhancer_db + "_prev" +
                        TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());

                StringBuilder sb_specific_db = new StringBuilder();
                sb_specific_db.append(
                                TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedFormat.get())
                        .append("\n");

                boolean needUplift = true;
                String ref_genome_before = "";

                try (BufferedReader dbReader = new BufferedReader(new FileReader(f_input_enhancerdbs_DB)))
                {
                    String line_db;
                    dbReader.readLine();

                    while ((line_db = dbReader.readLine()) != null)
                    {
                        String[] split = line_db.split("\t");
                        String species_refgenome = split[4];
                        String species_code = species_refgenome.replaceAll("[^A-Za-z]", "");

                        if (species_code.equals(igv_species_symbol))
                        {
                            if (species_refgenome.equals(TFPRIO.configs.igv.speciesReferenceGenome.get()))
                            {
                                sb_merged_dbs.append(line_db);
                                sb_merged_dbs.append("\n");
                                needUplift = false;
                                continue;
                            }

                            if (ref_genome_before.equals(""))
                            {
                                ref_genome_before = species_refgenome;
                            }

                            if (!species_refgenome.equals(ref_genome_before))
                            {
                                continue;
                            }

                            sb_specific_db.append(line_db);
                            sb_specific_db.append("\n");
                        }
                    }
                } catch (IOException e)
                {
                    e.printStackTrace();
                }

                if (needUplift)
                {
                    writeFile(f_output_filtered_db, sb_specific_db.toString());

                    File f_not_upliftedDB = extend(f_output_enhancer_dbs, enhancer_db + "_prev" +
                            TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());

                    File f_upliftedDB = new File(
                            f_output_enhancer_dbs.getAbsolutePath() + File.separator + enhancer_db +
                                    TFPRIO.configs.deSeq2.fileStructure.s_preprocessing_genePositions_mergedEnhancerDbs_bedSuffix.get());

                    File f_uplift_script_enhancer = extend(f_output_enhancer_dbs, enhancer_db + "_" +
                            TFPRIO.configs.deSeq2.fileStructure.f_preprocessing_genePositions_uplift.get().getName());

                    StringBuilder sb_pyuplift_db = new StringBuilder();

                    sb_pyuplift_db.append("import pip\n" + "\n" + "def import_or_install(package):\n" + "    try:\n" +
                            "        __import__(package)\n" + "    except ImportError:\n" +
                            "        pip.main(['install', package])\n" + "\n" + "import_or_install(\"pyliftover\")\n" +
                            "import_or_install(\"pandas\")\n" + "import_or_install(\"re\")\n" + "\n" + "import re\n" +
                            "import pandas as pd\n" + "from pyliftover import LiftOver\n" + "\n" +
                            "def convert(c, x, y, converter):\n" +
                            "    first = converter.convert_coordinate(c, int(x))\n" +
                            "    second = converter.convert_coordinate(c, int(y))\n" + "\n" +
                            "    if(first is None or second is None):\n" + "        return None, None\n" + "\n" +
                            "    if len(first) == 0 or len(second) == 0:\n" + "        return None, None\n" + "\n" +
                            "    return str(first[0][1]), str(second[0][1])\n" + "\n" + "def main():\n" +
                            "    path_to_X = \"" + f_output_filtered_db.getAbsolutePath() + "\"\n" +
                            "    path_to_newSave = \"" + f_upliftedDB.getAbsolutePath() + "\"\n" + "\n" +
                            "    version_of_our_dat=\"" + TFPRIO.configs.igv.speciesReferenceGenome.get() + "\"\n" +
                            "    version_of_their_dat" + "=\"" + ref_genome_before + "\"\n" + "\n" +
                            "    please_convert = True\n" + "    if(version_of_our_dat == version_of_their_dat):\n" +
                            "        please_convert=False\n" + "\n" + "\n" +
                            "    df = pd.read_csv(path_to_X, sep=\"\\t\")\n" + "\n" + "    data_output = []\n" + "\n" +
                            "    column_names=[]\n" + "\n" + "    for col in df.columns:\n" +
                            "        column_names.append(col)\n" + "\n" +
                            "    converter = LiftOver(version_of_their_dat, version_of_our_dat)\n" + "\n" +
                            "    for index, row in df.iterrows():\n" +
                            "        mgi_symbol = row[column_names.__getitem__(3)]\n" +
                            "        #chromosome_not_edited=str(row[column_names.__getitem__(0)])\n" +
                            "        chromosome = str(row[column_names.__getitem__(0)])\n" +
                            "        start_position = row[column_names.__getitem__(1)]\n" +
                            "        end_position = row[column_names.__getitem__(2)]\n" + "\n" +
                            "        ref_genome = row[column_names.__getitem__(4)]\n" + "\n" +
                            "        search_string = chromosome+\":\"+str(start_position)+\"-\"+str(end_position)\n" +
                            "        x_converted=start_position\n" + "        y_converted=end_position\n" + "\n" +
                            "        if please_convert:\n" +
                            "            x_converted, y_converted = convert(chromosome, start_position, end_position, converter)\n" +
                            "            if(x_converted is None or y_converted is None):\n" +
                            "                continue\n" + "\n" + "        start_position=x_converted\n" +
                            "        end_position=y_converted\n" + "\n" +
                            "        data_output.append([chromosome,start_position,end_position,mgi_symbol,'" +
                            TFPRIO.configs.igv.speciesReferenceGenome.get() + "'])\n" + "\n" +
                            "        #string_converted = chromosome+\":\"+x_converted+\"-\"+y_converted\n" +
                            "        #print(\"X\")\n" + "\n" +
                            "    df_converted = pd.DataFrame(data_output, columns=column_names)\n" +
                            "    df_converted.to_csv(path_to_newSave,sep=\"\\t\",index=False)\n" + "\n" + "main()\n");


                    BufferedWriter bw_pyuplift_db = new BufferedWriter(new FileWriter(f_uplift_script_enhancer));
                    bw_pyuplift_db.write(sb_pyuplift_db.toString());
                    bw_pyuplift_db.close();

                    logger.info("Uplift db " + enhancer_db + " positions to correct genome version");

                    String command_pyuplift_db = "python3 " + f_uplift_script_enhancer.getAbsolutePath();

                    logger.info("executing command: " + command_pyuplift_db);

                    Process child_pyuplift_db = Runtime.getRuntime().exec(command_pyuplift_db);
                    int code_pyuplift_db = child_pyuplift_db.waitFor();
                    switch (code_pyuplift_db)
                    {
                        case 0:
                            break;
                        case 1:
                            String message = child_pyuplift_db.getErrorStream().toString();
                            throw new Exception(message);
                    }

                    BufferedReader br_input_uplifted_names =
                            new BufferedReader(new FileReader(f_upliftedDB.getAbsolutePath()));
                    String line_input_uplifted_names = br_input_uplifted_names.readLine();
                    while ((line_input_uplifted_names = br_input_uplifted_names.readLine()) != null)
                    {
                        sb_merged_dbs.append(line_input_uplifted_names);
                        sb_merged_dbs.append("\n");
                    }
                    br_input_uplifted_names.close();
                }
            }

            File f_dbs_merged = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                    options_intern.file_suffix_deseq2_preprocessing_gene_positions_merged_enhancer_dbs);

            BufferedWriter bw_merged_enhancer_dbs = new BufferedWriter(new FileWriter(f_dbs_merged));
            bw_merged_enhancer_dbs.append(sb_merged_dbs.toString());
            bw_merged_enhancer_dbs.close();
        }
         */

        logger.info("Finished creating gene positions for all RNA-seq data.");
    }
}
