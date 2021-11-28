package tcga_atac_seq_helpers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

class download_building_blocks
{
    public String input_directory_mapping_sheet = "";
    String name = "";

}

public class Check_downloads
{

    public static void main(String[] args) throws IOException
    {
        download_building_blocks blocks_brca = new download_building_blocks();
        {
            blocks_brca.input_directory_mapping_sheet =
                    "/home/markus/data/COM2POSE/08_TCGA_ATAC_SEQ/03_BRCA/ATAC_SEQ/mapping_samples.csv";
            blocks_brca.name = "03_BRCA";

            check_double_sample_names(blocks_brca);
        }
    }

    private static void check_double_sample_names(download_building_blocks blocks_brca) throws IOException
    {
        System.out.println("CHECKING " + blocks_brca.name);

        File f_input_mapping_sheet = new File(blocks_brca.input_directory_mapping_sheet);
        BufferedReader br_mapping_sheet = new BufferedReader(new FileReader(f_input_mapping_sheet));
        String line_mapping_sheet = br_mapping_sheet.readLine();
        String[] header_mapping_sheet = line_mapping_sheet.split("\t");

        HashSet<String> sample_names = new HashSet<>();

        while ((line_mapping_sheet = br_mapping_sheet.readLine()) != null)
        {
            String[] split = line_mapping_sheet.split("\t");
            String sample_name = split[1];
            if (!sample_names.contains(sample_name))
            {
                sample_names.add(sample_name);
            } else
            {
                System.out.println("FOUND TWICE: " + sample_name);
            }
        }

        br_mapping_sheet.close();


    }

}
