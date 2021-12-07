package conversion_helpers;

import org.apache.commons.cli.*;
import util.Options_intern;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class RNA_seq_conversion
{

    public static void main(String[] args) throws Exception
    {
        Options_intern options_intern = new Options_intern();

        parseArguments(args, options_intern);

        if(options_intern.rna_seq_input_format.equals("ALL_IN_ONE"))
        {
            convert_all_together_to_seperate_gene_symbols(options_intern);
        }
        if(options_intern.rna_seq_input_format.equals("SALMON"))
        {
            convert_to_TF_PRIO_input(options_intern);
        }
    }

    private static void convert_to_TF_PRIO_input(Options_intern options_intern) throws IOException
    {
        File f_input = new File(options_intern.rna_seq_conversion_input);
        File f_output_root = new File(f_input.getParentFile().getAbsolutePath());

        //read in file
        ArrayList<String> ENSG_order = new ArrayList<>();
        ArrayList<String> sample_order = new ArrayList<>();
        HashMap<String,HashMap<String,Double>> sample_ensg_count= new HashMap<>();

        BufferedReader br_salmon_file =
                new BufferedReader(new FileReader(f_input));
        String line_salmon_file = br_salmon_file.readLine();
        String[] split_header_salmon_file = line_salmon_file.split("\t");
        for(int i = 2; i < split_header_salmon_file.length; i++)
        {
            HashMap<String,Double> sample_group = new HashMap<>();
            sample_ensg_count.put(split_header_salmon_file[i],sample_group);

            sample_order.add(split_header_salmon_file[i]);
        }
        while((line_salmon_file= br_salmon_file.readLine())!=null)
        {
            String[] split = line_salmon_file.split("\t");

            ENSG_order.add(split[0]);

            int position_in_sample_order_diff = -2;

            for(int i = 2; i < split.length; i++)
            {
                int position_in_sample_order = i + position_in_sample_order_diff;
                String sample_name = sample_order.get(position_in_sample_order);

                HashMap<String,Double> sample_hash = sample_ensg_count.get(sample_name);
                sample_hash.put(split[0],Double.parseDouble(split[i]));
            }
        }

        br_salmon_file.close();

        //create output directory
        File f_output_out_root = new File(f_output_root.getAbsolutePath()+File.separator+"formatted_input");
        f_output_out_root.mkdir();

        BufferedWriter bw_gene_list =
                new BufferedWriter(new FileWriter(new File(f_output_out_root.getAbsolutePath()+File.separator+
                        "Gene_id.txt")));
        bw_gene_list.write("Geneid\n");
        for(String key_ensg:ENSG_order)
        {
            bw_gene_list.write(key_ensg+"\n");
        }
        bw_gene_list.close();

        for(String key_group : options_intern.rna_seq_groups)
        {
            File f_out_group = new File(f_output_out_root.getAbsolutePath()+File.separator+key_group);
            f_out_group.mkdir();

            for(String key_sample : sample_ensg_count.keySet())
            {
                if(key_sample.matches(".*"+key_group+".*"))
                {
                    HashMap<String,Double> sample_to_counts = sample_ensg_count.get(key_sample);

                    BufferedWriter bw_group_counts =
                            new BufferedWriter(new FileWriter(f_out_group.getAbsolutePath()+File.separator+key_sample+".txt"));
                    bw_group_counts.write(key_sample+"\n");

                    for(String key_ensg: ENSG_order)
                    {
                        double count = sample_to_counts.get(key_ensg);
                        int count_i = 0;

                        if(options_intern.rna_seq_cut_to_integer)
                        {
                            count = Math.round(count);
                            count_i = (int) count;
                            bw_group_counts.write(count_i+"\n");
                        }
                        else
                        {
                            bw_group_counts.write(count+"\n");
                        }

                    }
                    bw_group_counts.close();
                }
            }



        }

        System.out.println("X");

    }

    private static void convert_all_together_to_seperate_gene_symbols(Options_intern options_intern) throws Exception
    {

        ArrayList<String> gene_names = new ArrayList<>();
        ArrayList<Boolean> gene_names_should_write = new ArrayList<>();

        //1. get gene_names order of one file to check for ENSG mapping
        File f_input_dir = new File(options_intern.rna_seq_conversion_input);
        File f_out = new File(options_intern.rna_seq_conversion_input + File.separator + "formatted_names.txt");


        for (File rna_seq_files : f_input_dir.listFiles())
        {
            if (rna_seq_files.isDirectory() || rna_seq_files.getName().matches(".*formatted.*") ||
                    rna_seq_files.getName().matches(".*gene_mapping.*"))
            {
                continue;
            }

            BufferedReader br = new BufferedReader(new FileReader(rna_seq_files));
            String line = "";
            int count_genes = 0;
            while ((line = br.readLine()) != null)
            {
                if (line.startsWith("#") || line.startsWith("_"))
                {
                    continue;
                }

                String[] split = line.split("\t");


                gene_names.add(split[0]);
                count_genes++;
            }
            br.close();
            break;

        }

        //count files to be processed
        int count_files_to_be_processed = 0;
        for (File rna_seq_files : f_input_dir.listFiles())
        {
            if (rna_seq_files.isDirectory() || rna_seq_files.getName().matches(".*formatted.*") ||
                    rna_seq_files.getName().matches(".*gene_mapping.*"))
            {
                continue;
            }
            count_files_to_be_processed++;
        }

        HashMap<String, String> symbol_ensg = new HashMap<>();

        //perform ENSG mapping
        if (!options_intern.rna_seq_conversion_biomart_species_name.equals(""))
        {
            //3. if biomart option is set:
            //a) convert symbol to ENSG

            File f_result_script = new File(f_input_dir.getAbsolutePath() + File.separator + "gene_mapping.R");
            File f_result_r = new File(f_input_dir.getAbsolutePath() + File.separator + "gene_mapping.csv");
            StringBuilder sb = new StringBuilder();
            sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                    "  install.packages(\"BiocManager\")\n" + "\n" +
                    "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" + "  BiocManager::install(\"biomaRt\")\n" +
                    "\n" + "library('biomaRt')\n");

            /*
            mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror="asia")
            df$id <- NA
            G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=df$Geneid,mart= mart)
             */

            sb.append("df <- read.csv('" + f_out.getAbsolutePath() + "')\n");

            sb.append("mart <- useEnsembl(biomart=\"ENSEMBL_MART_ENSEMBL\",dataset=\"" +
                    options_intern.rna_seq_conversion_biomart_species_name + "\", mirror=\"asia\")");
            sb.append("\n");

            sb.append(
                    "df$id <- NA\n" + "G_list <- getBM(filters= \"" + options_intern.rna_seq_conversion_biomart_column +
                            "\", attributes= c(\"ensembl_gene_id\",\"" +
                            options_intern.rna_seq_conversion_biomart_column + "\"),values=df$Geneid,mart= mart)\n" +
                            "write.table(G_list,\"" + f_result_r.getAbsolutePath() +
                            "\", row.names = FALSE, quote = F, sep=\"\\t\")\n");

            BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_result_script.getAbsolutePath()));
            bw_script.write(sb.toString());
            bw_script.close();

            //b) execute R script
/*
            String command = "Rscript " + f_result_script.getAbsolutePath();
            Process child = Runtime.getRuntime().exec(command);
            System.out.println("[R-SCRIPT] Running script " + f_result_script.getName()+": " + command);
            int code = child.waitFor();
            switch (code) {
                case 0:
                    break;
                case 1:
                    String message = child.getErrorStream().toString();
                    throw new Exception(message);
            }*/

            //c) write new gene_file gene_names file
            BufferedReader br_mapping = new BufferedReader(new FileReader(f_result_r));
            String line_mapping = br_mapping.readLine();
            while ((line_mapping = br_mapping.readLine()) != null)
            {
                String[] split = line_mapping.split("\t");
                symbol_ensg.put(split[1], split[0]);
            }

            br_mapping.close();

            HashSet<String> already_founds_ensgs = new HashSet<>();


            BufferedWriter bw_formatted = new BufferedWriter(new FileWriter(f_out));
            bw_formatted.write("Geneid");
            bw_formatted.newLine();
            int count = 0;
            for (String s : gene_names)
            {
                if (!symbol_ensg.containsKey(s))
                {
                    gene_names_should_write.add(false);
                } else
                {

                    if (already_founds_ensgs.contains(s))
                    {
                        gene_names_should_write.add(false);
                    } else
                    {
                        bw_formatted.write(symbol_ensg.get(s));
                        bw_formatted.newLine();
                        gene_names_should_write.add(true);
                        already_founds_ensgs.add(s);
                    }
                }
                count++;
            }
            bw_formatted.close();

            /*
            int row_count=0;
            for(Boolean b: gene_names_should_write)
            {
                if(b)
                    row_count++;
            }

            HashMap<String,Integer> ensg_position = new HashMap<>();
            ArrayList<String> mapped_list = new ArrayList<>();
            ArrayList<String> is_double_or_not_found = new ArrayList<>();
            HashSet<String> already_found_ensgs_count = new HashSet<>();
            int count_doubles = 0;
            int count_position=0;
            for(String names : gene_names)
            {
                if(symbol_ensg.containsKey(names))
                {
                    String name_mapped = symbol_ensg.get(names);
                    if(already_found_ensgs_count.contains(name_mapped))
                    {
                        is_double_or_not_found.add("D");
                        count_doubles++;
                    }
                    else
                    {
                        is_double_or_not_found.add("T");
                        already_found_ensgs_count.add(name_mapped);
                    }
                }
                else
                {
                    is_double_or_not_found.add("F");
                    continue;
                }

                count_position++;

            }

            String[][] files_to_counts = new String[row_count-count_doubles][count_files_to_be_processed+1];*/


        } else
        {
            BufferedWriter bw_ensg_out = new BufferedWriter(new FileWriter(f_out));
            bw_ensg_out.write("Geneid");
            bw_ensg_out.newLine();
            for (String gene_name : gene_names)
            {
                bw_ensg_out.write(gene_name);
                bw_ensg_out.newLine();
            }
            bw_ensg_out.close();

            for (String s : gene_names)
            {
                gene_names_should_write.add(true);
            }


        }

        //2. write gene counts in right order
        int count_files = 0;
        for (File rna_seq_files : f_input_dir.listFiles())
        {
            if (rna_seq_files.isDirectory() || rna_seq_files.getName().matches(".*formatted.*") ||
                    rna_seq_files.getName().matches(".*gene_mapping.*"))
            {
                continue;
            }

            File f_output = new File(rna_seq_files.getAbsolutePath() + "_formatted.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(f_output));
            bw.write(rna_seq_files.getName());
            bw.newLine();
            BufferedReader br = new BufferedReader(new FileReader(rna_seq_files));
            String line = "";
            int count_genes = 0;
            while ((line = br.readLine()) != null)
            {
                if (line.startsWith("#") || line.startsWith("_"))
                {
                    continue;
                }

                String[] split = line.split("\t");
                if (gene_names_should_write.get(count_genes))
                {
                    bw.write(split[1]);
                    bw.newLine();
                }
                count_genes++;
            }
            br.close();
            bw.close();

            count_files++;

        }
    }

    private static void parseArguments(String[] args, Options_intern options_intern)
    {
        Options options = new Options();

        Option opt_working_dir =
                new Option("i", "file-directory", true, "[REQ]: file-directory with RNA-seq data in other format");
        opt_working_dir.setRequired(true);
        options.addOption(opt_working_dir);

        Option opt_input_format = new Option("f","format",true,"[REQ]: ALL_IN_ONE if first column in gene and the " +
                "rest is just one sample after another, SALMON if you provide salmon.merged.gene_counts.tsv from " +
                "nfcore or salmon");
        opt_input_format.setRequired(true);
        options.addOption(opt_input_format);

        Option opt_groups = new Option("g","groups",true,"[REQ for SALMON format]: write groups like that: untreated;" +
                "IFNb");
        options.addOption(opt_groups);

        Option opt_cut_to_integer = new Option("I","cut-to-integer",false,"[OPT]: gene counts are cut on their " +
                "decimal to integers. Needed for DESeq2.");
        options.addOption(opt_cut_to_integer);

        Option opt_biomart_species = new Option("s", "biomart-species", true,
                "[OPT]: if biomart-species is set, automatic transfer to ENSG will be performed, NOTE: must be set if gene names are in SYMBOL form");
        //opt_biomart_species.setRequired(true);
        options.addOption(opt_biomart_species);
        Option opt_biomart_column = new Option("c", "biomart-column", true,
                "[OPT]: if biomart-species is set, the column name must be set");
        //opt_biomart_species.setRequired(true);
        options.addOption(opt_biomart_column);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try
        {
            cmd = parser.parse(options, args);


            options_intern.rna_seq_conversion_input = cmd.getOptionValue("file-directory");
            options_intern.rna_seq_input_format = cmd.getOptionValue("format");
            options_intern.rna_seq_conversion_biomart_species_name = cmd.getOptionValue("biomart-species", "");
            options_intern.rna_seq_conversion_biomart_column = cmd.getOptionValue("biomart-column", "");

            if (cmd.hasOption("cut-to-integer"))
            {
                options_intern.rna_seq_cut_to_integer = true;
            }

            if (cmd.hasOption("groups"))
            {
                String groups = cmd.getOptionValue("groups");
                String[] split_groups = groups.split(";");
                Collections.addAll(options_intern.rna_seq_groups, split_groups);
            }


        } catch (ParseException e)
        {
            System.out.println(e.getMessage());
            formatter.printHelp("-i <file-directory> -f <format> [-I] [-s <biomart-species>] [-c <biomart-column>]",
                    options);
            System.exit(1);
        }


    }
}
