package tcga_atac_seq_helpers;

import java.io.*;

public class Create_environments
{

    public static void main(String[] args) throws IOException
    {
        String root_file_input = "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ";
        String root_file_working_directories = "/home/markus/data/COM2POSE/01C_TCGA";
        String default_conifg =
                "/home/markus/data/COM2POSE/01C_TCGA/01_TGCT/default_run_inc_tgene/com2pose_template.cfg";
        String default_command = "/home/markus/data/COM2POSE/01C_TCGA/01_TGCT/default_run_inc_tgene/commandline.txt";

        create_and_copy_config_and_commandline(root_file_input, default_command, default_conifg,
                root_file_working_directories);
    }

    private static void create_and_copy_config_and_commandline(String root_file, String default_command,
                                                               String default_conifg,
                                                               String root_file_working_directories) throws IOException
    {

        File f_root_file_working_directories = new File(root_file_working_directories);
        File f_root = new File(root_file);

        File f_default_config = new File(default_conifg);
        File f_default_command = new File(default_command);

        for (File f_file : f_root.listFiles())
        {
            if (f_file.isDirectory())
            {
                String name = f_file.getName();

                if (name.startsWith("00_requiredData"))
                {
                    continue;
                }

                //check if two or more subtypes are available
                //if yes -> continue;
                //if no -> break;

                String atac_seq = "ATAC_SEQ";
                String sorted_data = "02_sorted_data";

                File f_dir_to_atac_seq_data = new File(
                        f_file.getAbsoluteFile() + File.separator + atac_seq + File.separator + atac_seq +
                                File.separator + sorted_data);
                if (!f_dir_to_atac_seq_data.exists())
                {
                    System.out.println("[ERROR]: ATAC_SEQ file does not exist");
                    System.out.println("[ERROR]: at " + f_file.getName());
                    continue;
                    //System.exit(1);
                }
                String rna_seq = "RNA_SEQ";
                String formatted = "03_formatted";

                File f_dir_rna_seq_data =
                        new File(f_file.getAbsolutePath() + File.separator + rna_seq + File.separator + formatted);
                if (!f_dir_rna_seq_data.exists())
                {
                    System.out.println("[ERROR]: RNA_SEQ file does not exist");
                    System.out.println("[ERROR]: at " + f_file.getName());
                    continue;
                    //System.exit(1);
                }

                File f_dir_rna_seq_data_geneids =
                        new File(f_dir_rna_seq_data.getAbsolutePath() + File.separator + "names_formatted.txt");
                if (!f_dir_rna_seq_data_geneids.exists())
                {
                    System.out.println("[ERROR]: RNA_SEQ geneids file does not exist");
                    System.exit(1);
                }

                if (f_dir_to_atac_seq_data.listFiles().length > 1)
                {
                    //create head folder
                    File f_head = new File(f_root_file_working_directories.getAbsolutePath() + File.separator + name);
                    f_head.mkdir();

                    String name_default = "default_run";
                    String name_default_inc_tgene = "default_run_inc_tgene";

                    //create two subfolder -> inc tgene and no tgene
                    File f_default = new File(f_head.getAbsolutePath() + File.separator + name_default);
                    f_default.mkdir();

                    File f_default_inc_tgene =
                            new File(f_head.getAbsolutePath() + File.separator + name_default_inc_tgene);
                    f_default_inc_tgene.mkdir();

                    //in each create working_dir folder
                    File f_default_working_dir = new File(f_default.getAbsolutePath() + File.separator + "working_dir");
                    f_default_working_dir.mkdir();

                    File f_default_inc_tgene_working_dir =
                            new File(f_default_inc_tgene.getAbsolutePath() + File.separator + "working_dir");
                    f_default_inc_tgene_working_dir.mkdir();

                    File f_default_config_file =
                            new File(f_default.getAbsolutePath() + File.separator + f_default_config.getName());
                    File f_default_config_file_inc_tgene = new File(
                            f_default_inc_tgene.getAbsolutePath() + File.separator + f_default_config.getName());

                    //create commandline.txt

                    StringBuilder sb_default_config = new StringBuilder();
                    sb_default_config.append("-c\n");
                    sb_default_config.append(f_default_config_file.getAbsolutePath());
                    sb_default_config.append("\n-w\n");
                    sb_default_config.append(f_default_working_dir.getAbsolutePath());
                    sb_default_config.append("\n-p\n");
                    sb_default_config.append("/home/markus/IdeaProjects/COM2POSE");

                    StringBuilder sb_default_config_inc_tgene = new StringBuilder();
                    sb_default_config_inc_tgene.append("-c\n");
                    sb_default_config_inc_tgene.append(f_default_config_file_inc_tgene.getAbsolutePath());
                    sb_default_config_inc_tgene.append("\n-w\n");
                    sb_default_config_inc_tgene.append(f_default_inc_tgene_working_dir.getAbsolutePath());
                    sb_default_config_inc_tgene.append("\n-p\n");
                    sb_default_config_inc_tgene.append("/home/markus/IdeaProjects/COM2POSE");
                    sb_default_config_inc_tgene.append("\n-t\n");
                    sb_default_config_inc_tgene.append("/home/markus/meme");

                    File f_default_command_txt =
                            new File(f_default.getAbsolutePath() + File.separator + "commandline.txt");
                    BufferedWriter bw_default_command = new BufferedWriter(new FileWriter(f_default_command_txt));
                    bw_default_command.write(sb_default_config.toString());
                    bw_default_command.close();

                    File f_default_command_txt_inctgene =
                            new File(f_default_inc_tgene.getAbsolutePath() + File.separator + "commandline.txt");
                    BufferedWriter bw_default_command_inctgene =
                            new BufferedWriter(new FileWriter(f_default_command_txt_inctgene));
                    bw_default_command_inctgene.write(sb_default_config_inc_tgene.toString());
                    bw_default_command_inctgene.close();

                    //create config with current paths
                    StringBuilder sb_altered_config_file = new StringBuilder();
                    BufferedReader br_old_config_file = new BufferedReader(new FileReader(f_default_config));
                    String line_old_config_file = "";
                    while ((line_old_config_file = br_old_config_file.readLine()) != null)
                    {
                        if (line_old_config_file.startsWith("deseq2_input_directory"))
                        {
                            sb_altered_config_file.append("deseq2_input_directory=\"");
                            sb_altered_config_file.append(f_dir_rna_seq_data.getAbsolutePath());
                            sb_altered_config_file.append("\"\n");

                        } else if (line_old_config_file.startsWith("deseq2_input_gene_id"))
                        {
                            sb_altered_config_file.append("deseq2_input_gene_id=\"");
                            sb_altered_config_file.append(f_dir_rna_seq_data_geneids.getAbsolutePath());
                            sb_altered_config_file.append("\"\n");

                        } else if (line_old_config_file.startsWith("tepic_input_directory"))
                        {
                            sb_altered_config_file.append("tepic_input_directory=\"");
                            sb_altered_config_file.append(f_dir_to_atac_seq_data.getAbsolutePath());
                            sb_altered_config_file.append("\"\n");

                        } else
                        {
                            sb_altered_config_file.append(line_old_config_file);
                            sb_altered_config_file.append("\n");
                        }

                    }
                    br_old_config_file.close();

                    //!!DO NOT OVERWRITE MODIFIED TEMPLATES!!//
                    /*
                    BufferedWriter bw_config_inc_tgene = new BufferedWriter(new FileWriter(f_default_config_file_inc_tgene));
                    bw_config_inc_tgene.write(sb_altered_config_file.toString());
                    bw_config_inc_tgene.close();

                    BufferedWriter bw_config= new BufferedWriter(new FileWriter(f_default_config_file));
                    bw_config.write(sb_altered_config_file.toString());
                    bw_config.close();*/


                }
            }
        }
        System.out.println("X");
    }
}
