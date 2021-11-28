package tcga_atac_seq_helpers;

import java.util.HashSet;

public class atac_seq_download_cancer_blocks
{

    public String input_gdc_client = "/home/markus/Downloads/gdc-client_v1.6.1_Ubuntu_x64/gdc-client";
    public String input_token = "/home/markus/Downloads/gdc-user-token.2021-07-16T10_22_43.173Z.txt";

    // these ones always change!
    public String name = "";

    public String input_manifest = "";
    public String input_gdc_sample_sheet = "";
    public String input_clinical = "";

    public String output_directory = "";

    public HashSet<String> groups = new HashSet<>();
    public String other_samples = "";

    //rna_seq only
    public String input_dir_root = "";
    public int number_groups = -1;


}
