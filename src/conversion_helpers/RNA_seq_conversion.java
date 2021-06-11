package conversion_helpers;

import org.apache.commons.cli.*;
import util.Options_intern;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class RNA_seq_conversion {

    public static void main(String[] args) throws Exception {
        Options_intern options_intern= new Options_intern();

        parseArguments(args, options_intern);

        convert_all_together_to_seperate_gene_symbols(options_intern);
    }

    private static void convert_all_together_to_seperate_gene_symbols(Options_intern options_intern) throws Exception {

        ArrayList<String> gene_names = new ArrayList<>();
        ArrayList<Boolean> gene_names_should_write = new ArrayList<>();

        //1. get gene_names order of one file to check for ENSG mapping
        File f_input_dir = new File(options_intern.rna_seq_conversion_input);
        File f_out = new File(options_intern.rna_seq_conversion_input+File.separator+"formatted_names.txt");


        for(File rna_seq_files : f_input_dir.listFiles())
        {
            if (rna_seq_files.isDirectory() || rna_seq_files.getName().matches(".*formatted.*") || rna_seq_files.getName().matches(".*gene_mapping.*")) {
                continue;
            }

            BufferedReader br = new BufferedReader(new FileReader(rna_seq_files));
            String line ="";
            int count_genes = 0;
            while((line = br.readLine())!=null)
            {
                if(line.startsWith("#")||line.startsWith("_"))
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

        //perform ENSG mapping
        if(!options_intern.rna_seq_conversion_biomart_species_name.equals("")) {
            //3. if biomart option is set:
            //a) convert symbol to ENSG

            File f_result_script = new File(f_input_dir.getAbsolutePath() + File.separator + "gene_mapping.R");
            File f_result_r = new File(f_input_dir.getAbsolutePath() + File.separator + "gene_mapping.csv");
            StringBuilder sb = new StringBuilder();
            sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                    "  install.packages(\"BiocManager\")\n" +
                    "\n" +
                    "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" +
                    "  BiocManager::install(\"biomaRt\")\n" +
                    "\n" +
                    "library('biomaRt')\n");

            /*
            mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror="asia")
            df$id <- NA
            G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=df$Geneid,mart= mart)
             */

            sb.append("df <- read.csv('" + f_out.getAbsolutePath() + "')\n");

            sb.append("mart <- useEnsembl(biomart=\"ENSEMBL_MART_ENSEMBL\",dataset=\"" + options_intern.rna_seq_conversion_biomart_species_name + "\", mirror=\"asia\")");
            sb.append("\n");

            sb.append("df$id <- NA\n" +
                    "G_list <- getBM(filters= \"" + options_intern.rna_seq_conversion_biomart_column + "\", attributes= c(\"ensembl_gene_id\",\"" + options_intern.rna_seq_conversion_biomart_column + "\"),values=df$Geneid,mart= mart)\n" +
                    "write.table(G_list,\"" + f_result_r.getAbsolutePath() + "\", row.names = FALSE, quote = F, sep=\"\\t\")\n");

            BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_result_script.getAbsolutePath()));
            bw_script.write(sb.toString());
            bw_script.close();

            //b) execute R script

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
            }

            //c) write new gene_file gene_names file
            HashMap<String, String> symbol_ensg = new HashMap<>();
            BufferedReader br_mapping = new BufferedReader(new FileReader(f_result_r));
            String line_mapping = br_mapping.readLine();
            while ((line_mapping = br_mapping.readLine()) != null) {
                String[] split = line_mapping.split("\t");
                symbol_ensg.put(split[1], split[0]);
            }

            br_mapping.close();


            BufferedWriter bw_formatted = new BufferedWriter(new FileWriter(f_out));
            bw_formatted.write("Geneid");
            bw_formatted.newLine();
            for (String s : gene_names) {
                if (!symbol_ensg.containsKey(s)) {
                    gene_names_should_write.add(false);
                }
                else {
                    bw_formatted.write(symbol_ensg.get(s));
                    bw_formatted.newLine();
                    gene_names_should_write.add(true);
                }
            }
            bw_formatted.close();
        }
        else
        {
            BufferedWriter bw_ensg_out = new BufferedWriter(new FileWriter(f_out));
            bw_ensg_out.write("Geneid");
            bw_ensg_out.newLine();
            for(String gene_name : gene_names)
            {
                bw_ensg_out.write(gene_name);
                bw_ensg_out.newLine();
            }
            bw_ensg_out.close();

            for(String s: gene_names)
            {
                gene_names_should_write.add(true);
            }
        }

        //2. write gene counts in right order
        int count_files = 0;
        for(File rna_seq_files : f_input_dir.listFiles())
        {
            if(rna_seq_files.isDirectory() || rna_seq_files.getName().matches(".*formatted.*") || rna_seq_files.getName().matches(".*gene_mapping.*"))
            {
                continue;
            }

            File f_output = new File(rna_seq_files.getAbsolutePath()+"_formatted.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(f_output));
            bw.write(rna_seq_files.getName());
            bw.newLine();
            BufferedReader br = new BufferedReader(new FileReader(rna_seq_files));
            String line ="";
            int count_genes = 0;
            while((line = br.readLine())!=null)
            {
                if(line.startsWith("#")||line.startsWith("_"))
                {
                    continue;
                }

                String[] split = line.split("\t");
                if(gene_names_should_write.get(count_genes))
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




        System.out.println("X");


    }

    private static void parseArguments(String[] args, Options_intern options_intern)
    {
        Options options = new Options();

        Option opt_working_dir = new Option("i","file-directory",true,"[REQ]: file-directory with RNA-seq data in other format");
        opt_working_dir.setRequired(true);
        options.addOption(opt_working_dir);

        Option opt_biomart_species = new Option("s","biomart-species",true,"[OPT]: if biomart-species is set, automatic transfer to ENSG will be performed, NOTE: must be set if gene names are in SYMBOL form");
        //opt_biomart_species.setRequired(true);
        options.addOption(opt_biomart_species);
        Option opt_biomart_column = new Option("c","biomart-column",true,"[OPT]: if biomart-species is set, the column name must be set");
        //opt_biomart_species.setRequired(true);
        options.addOption(opt_biomart_column);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try
        {
            cmd = parser.parse(options, args);


            options_intern.rna_seq_conversion_input = cmd.getOptionValue("file-directory");
            options_intern.rna_seq_conversion_biomart_species_name=cmd.getOptionValue("biomart-species","");
            options_intern.rna_seq_conversion_biomart_column=cmd.getOptionValue("biomart-column","");



        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("-i <file-directory> [-s <biomart-species>] [-c <biomart-column>]",options);
            System.exit(1);
        }



    }
}
